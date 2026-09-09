function F_vec = computeFluxDerivative(U, fluxLimiter, lambda, gamma, dx, kappa, recon, CD_Term_Order, riemannSolver, xs, time, sourceTermFn)
% computeFluxDerivative - compute the spatial flux derivative for Euler equations.
%
% This function computes the right-hand side dU/dt for the conservative
% variables U using a selected reconstruction scheme and Riemann solver.
%
% Inputs:
%   U             - conservative variables (3 x N, including ghost cells)
%   fluxLimiter   - integer ID of the flux limiter (used only for MUSCL)
%   lambda        - MUSCL lambda scaling factor
%   gamma         - ratio of specific heats
%   dx            - grid spacing
%   kappa         - MUSCL parameter
%   recon         - reconstruction method ('ROE','MUSCL','WENO3','WENO5')
%   CD_Term_Order - order of the central differencing term
%   riemannSolver - Riemann solver name ('Roe','Roe-Pike','HLL','HLLC','Osher')
%   xs            - grid locations (3 x N-matching, required for some
%                   solvers and for sourceTermFn)
%   time          - current simulation time (required for some solvers)
%   sourceTermFn  - OPTIONAL function handle @(x) -> 3xM source term,
%                   evaluated at the interior grid locations and added to
%                   the RHS (for the method-of-manufactured-solutions
%                   verification path; the PDE solved there is
%                   dU/dt + dF/dx = S(x), i.e. dU/dt = S(x) - dF/dx).
%                   Omit (or pass []) for the standard, unforced Euler
%                   equations.
%
% Output:
%   F_vec - flux derivative (3 x N, including ghost cells)

if nargin < 12
    sourceTermFn = [];
end

% Select and apply reconstruction
switch recon
    case 'ROE'
        [URi,ULi,URim1,ULim1] = musclRec(U, 18, lambda, kappa);
    case 'MUSCL'
        [URi,ULi,URim1,ULim1] = musclRec(U, fluxLimiter, lambda, kappa);
    case 'WENO3'
        [URi,ULi] = weno3Rec(U(:,2:end));
        [URim1,ULim1] = weno3Rec(U(:,1:end-1));
    case 'WENO5'
        [URi,ULi] = weno5Rec(U(:,2:end));
        [URim1,ULim1] = weno5Rec(U(:,1:end-1));
    otherwise
        error('Unknown reconstruction method: %s', recon);
end

% Select and apply Riemann solver
switch riemannSolver
    case 'Roe'
        [Fip12, Fim12] = roeSolver(U, URi, ULi, URim1, ULim1, gamma, CD_Term_Order);
    case 'Roe-Pike'
        [Fip12, Fim12] = roePikeSolver(U, URi, ULi, URim1, ULim1, gamma, CD_Term_Order);
    case 'HLL'
        [Fip12, Fim12] = hllSolver(U, URi, ULi, URim1, ULim1, xs, gamma, time);
    case 'HLLC'
        [Fip12, Fim12] = hllcSolver(U, URi, ULi, URim1, ULim1, xs, gamma, time);
    case 'Osher'
        [Fip12, Fim12] = osherSolver(U, URi, ULi, URim1, ULim1, xs, gamma, time);
    otherwise
        error('Unknown Riemann solver: %s', riemannSolver);
end

% Compute flux derivative (including padding for ghost cells)
F_vec = -(1/dx) * (Fip12 - Fim12);
% Keep ghost cells unchanged (they are updated by boundary condition routines)
% FIXED: previously used the *total* column deficit (size(U,2)-size(F_vec,2))
% as the padding width on *each* side, doubling the total padding added
% (e.g. 6 ghost columns missing became 6+6=12 added instead of 3+3=6),
% so F_vec came out wider than U. Now splits the deficit evenly, matching
% the equal number of ghost cells on each side.
ghost = (size(U,2) - size(F_vec,2)) / 2;

if ~isempty(sourceTermFn)
    xint = xs(ghost+1:end-ghost);
    F_vec = F_vec + sourceTermFn(xint);
end

BCs = zeros(3, ghost);
F_vec = [BCs, F_vec, BCs];
end
