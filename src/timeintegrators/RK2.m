function U = RK2(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn)
% 2nd order Runge-Kutta
% sourceTermFn - OPTIONAL, see computeFluxDerivative.m

if nargin < 13
    sourceTermFn = [];
end

k1 = U + dt*computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U = 0.5*U + 0.5*k1 + 0.5*dt*computeFluxDerivative(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
end
