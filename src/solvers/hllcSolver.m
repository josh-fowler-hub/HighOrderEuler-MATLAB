function [Fip12,Fim12] = hllcSolver(U,URi,ULi,URim1,ULim1,xs,gamma,time)
% hllcSolver - HLLC approximate Riemann solver (Toro Ch. 10), vectorized
% over interfaces. Computed independently at both the i+1/2 interface
% (from ULi/URi) and the i-1/2 interface (from ULim1/URim1).
%
% FIXED: previously ignored URim1/ULim1 entirely and set Fim12 = Fip12
% (both computed from the current interface only), which zeroed out the
% finite-volume flux difference at every cell. Also previously used the
% plain two-wave HLL average formula in the "middle" region rather than
% genuine HLLC star-state fluxes (the contact-preserving "C" was absent).
%
% Inputs:
%   U      - full conservative state array (unused except for signature
%            uniformity with the other Riemann solvers)
%   URi    - right reconstructed state at interface i+1/2
%   ULi    - left reconstructed state at interface i+1/2
%   URim1  - right reconstructed state at interface i-1/2
%   ULim1  - left reconstructed state at interface i-1/2
%   xs      - grid locations (unused)
%   gamma   - ratio of specific heats
%   time    - current time (unused)
%
% Outputs:
%   Fip12, Fim12 - numerical fluxes on the right and left side of each cell

Fip12 = hllcFlux(ULi, URi, gamma);
Fim12 = hllcFlux(ULim1, URim1, gamma);
end

function F = hllcFlux(UL, UR, gamma)
% Genuine HLLC flux (Toro eq. 10.37-10.43) between left state UL and
% right state UR.

[rhoL, uL, pL, aL, ~, ~, ~, ~] = consToPrim(UL, gamma);
[rhoR, uR, pR, aR, ~, ~, ~, ~] = consToPrim(UR, gamma);

% Outer wave-speed bounds via the adaptive (Toro 9.3-9.5) estimate, same
% as used by the other approximate Riemann solvers in this codebase.
% FIXED: floored/ceilinged against the raw data wave speeds uL-aL/uR+aR
% (standard Toro/Einfeldt robust HLL(C) wave-speed estimate) - the
% adaptive estimate alone can under-bound the true wave fan at some
% interfaces, which destabilizes the scheme (see hllSolver.m for the
% numerical confirmation).
[ustar, astarL, astarR] = adaptiveRiemannSolver(rhoL, rhoR, uL, uR, aL, aR, pL, pR, gamma);
SL = min(uL - aL, ustar - astarL);
SR = max(uR + aR, ustar + astarR);

FL = Fvec(rhoL, uL, pL, gamma);
FR = Fvec(rhoR, uR, pR, gamma);

EL = UL(3,:); % total energy per unit volume
ER = UR(3,:);

% Contact (middle) wave speed, Toro eq. 10.37
SM = (pR - pL + rhoL.*uL.*(SL-uL) - rhoR.*uR.*(SR-uR)) ./ (rhoL.*(SL-uL) - rhoR.*(SR-uR));

% Star states, Toro eq. 10.39
coefL = rhoL.*(SL-uL)./(SL-SM);
UstarL = [coefL; coefL.*SM; coefL.*(EL./rhoL + (SM-uL).*(SM + pL./(rhoL.*(SL-uL))))];

coefR = rhoR.*(SR-uR)./(SR-SM);
UstarR = [coefR; coefR.*SM; coefR.*(ER./rhoR + (SM-uR).*(SM + pR./(rhoR.*(SR-uR))))];

FstarL = FL + SL.*(UstarL - UL);
FstarR = FR + SR.*(UstarR - UR);

F = zeros(size(UL));
maskL  = SL >= 0;
maskML = (SL < 0) & (SM >= 0);
maskMR = (SM < 0) & (SR >= 0);
maskR  = SR < 0;

F(:,maskL)  = FL(:,maskL);
F(:,maskML) = FstarL(:,maskML);
F(:,maskMR) = FstarR(:,maskMR);
F(:,maskR)  = FR(:,maskR);
end
