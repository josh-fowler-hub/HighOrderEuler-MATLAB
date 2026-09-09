function U = RK3(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn)
% 3rd order Runge-Kutta
% sourceTermFn - OPTIONAL, see computeFluxDerivative.m

if nargin < 13
    sourceTermFn = [];
end

k1 = U + dt*computeFluxDerivative(U,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
k2 = 0.75*U + 0.25*k1 + 0.25*dt*computeFluxDerivative(k1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U = (1/3)*U + (2/3)*k2 + (2/3)*dt*computeFluxDerivative(k2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
end
