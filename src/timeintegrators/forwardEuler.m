function U = forwardEuler(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn)
% forwardEuler - First-order explicit Euler time integration.
% Performs a single time step: U^{n+1} = U^n + dt * dU/dt.
%
% sourceTermFn - OPTIONAL, see computeFluxDerivative.m

if nargin < 13
    sourceTermFn = [];
end

% Compute the spatial derivative (flux residual)
F = computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);

% Advance the solution
U = U + dt * F;
end
