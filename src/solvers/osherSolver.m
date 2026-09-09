function [Fip12,Fim12] = osherSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time)
% osherSolver - Osher-style numerical flux for 1D Euler equations.
%
% This implementation follows an Osher-style flux construction using a
% Roe-linearization along a straight-line path in state space. It
% evaluates the dissipation term D = R * |Lambda| * R^{-1} * dU along
% that path (via 2-point Gauss quadrature) and uses it to compute the
% interface flux. Computed independently at both the i+1/2 interface
% (from ULi/URi) and the i-1/2 interface (from ULim1/URim1).
%
% FIXED: previously ignored URim1/ULim1 entirely and set Fim12 = Fip12
% (both computed from the current interface only), which zeroed out the
% finite-volume flux difference at every cell.
%
% VECTORIZED: previously built the full 3x3 matrix |A| at every path
% point via a per-interface for-loop (absEulerJacobian.m + a second
% per-interface loop for the matrix-vector product), which is by far the
% dominant cost of this solver under Octave (no scalar-loop JIT). This
% version instead evaluates D*dU directly via a closed-form wave-strength
% decomposition, decomposing the FULL conservative jump dU = UR-UL into
% the local eigenbasis at each path point:
%   dp_eff = (gamma-1)*(dU3 - u0*dU2 + 0.5*u0^2*dU1)
%   alpha1 = (dp_eff - a0*(dU2-u0*dU1)) / (2*a0^2)
%   alpha2 = dU1 - dp_eff/a0^2
%   alpha3 = (dp_eff + a0*(dU2-u0*dU1)) / (2*a0^2)
%   |A|*dU = alpha1*|lambda1|*K1 + alpha2*|lambda2|*K2 + alpha3*|lambda3|*K3
% This is an exact identity for ANY reference state (u0,a0) - not the
% Roe-average-specific shortcut used in roeSolver.m/roePikeSolver.m for
% their own (different) dU - verified numerically to match the original
% per-index matrix computation to floating-point precision (~1e-16)
% before being adopted here.
%
% Inputs (all 3xN matrices):
%   U     - full conservative state array (unused, kept for signature
%           uniformity with the other Riemann solvers)
%   URi   - right reconstructed state at interface i+1/2
%   ULi   - left reconstructed state at interface i+1/2
%   URim1 - right reconstructed state at interface i-1/2
%   ULim1 - left reconstructed state at interface i-1/2
%   xs, time - unused
%   gamma - ratio of specific heats
%
% Outputs:
%   Fip12, Fim12 - the numerical flux on the right and left side of the cell

Fip12 = osherFlux(ULi, URi, gamma);
Fim12 = osherFlux(ULim1, URim1, gamma);
end

function F = osherFlux(UL, UR, gamma)
% Convert to primitive variables
[rhoL, uL, pL, ~, ~, ~, ~, ~] = consToPrim(UL, gamma);
[rhoR, uR, pR, ~, ~, ~, ~, ~] = consToPrim(UR, gamma);

% Physical fluxes for left and right states
FL = Fvec(rhoL, uL, pL, gamma);
FR = Fvec(rhoR, uR, pR, gamma);

% Difference in conservative variables (fixed for both quadrature points)
dU = UR - UL;
d1 = dU(1,:); d2 = dU(2,:); d3 = dU(3,:);

ones_row = ones(1, size(UL,2));

% Two-point Gauss quadrature for the path integral. Each point
% contributes 0.25*|A(path point)|*dU to Fdiff (matching the original
% D = 0.5*(absA(s1)+absA(s2)); Fdiff = 0.5*D*dU scaling exactly).
s = [0.211324865405187, 0.788675134594813];
Fdiff = zeros(size(UL));
for si = s
    Upath = UL + si .* dU;
    [~, uP, ~, aP, HP, ~, ~, ~] = consToPrim(Upath, gamma);

    lam1 = uP - aP;
    lam2 = uP;
    lam3 = uP + aP;

    K1 = [ones_row; lam1; HP - uP.*aP];
    K2 = [ones_row; uP; 0.5.*uP.^2];
    K3 = [ones_row; lam3; HP + uP.*aP];

    dp_eff = (gamma-1).*(d3 - uP.*d2 + 0.5.*uP.^2.*d1);
    alpha1 = (dp_eff - aP.*(d2 - uP.*d1)) ./ (2.*aP.^2);
    alpha2 = d1 - dp_eff./aP.^2;
    alpha3 = (dp_eff + aP.*(d2 - uP.*d1)) ./ (2.*aP.^2);

    Fdiff = Fdiff + 0.25.*(alpha1.*abs(lam1).*K1 + alpha2.*abs(lam2).*K2 + alpha3.*abs(lam3).*K3);
end

% Final flux
F = 0.5 .* (FL + FR) - Fdiff;
end
