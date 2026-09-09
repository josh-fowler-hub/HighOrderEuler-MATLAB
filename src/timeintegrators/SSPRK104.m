function U = SSPRK104(U,flux,lam,gamma,dx,dt,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn)
% SSPRK(10,4) strong stability preserving Runge-Kutta
% sourceTermFn - OPTIONAL, see computeFluxDerivative.m

if nargin < 13
    sourceTermFn = [];
end

U0 = U;
U1 = U0 + dt/6 *computeFluxDerivative(U0, flux, lam, gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U2 = U1 + dt/6 * computeFluxDerivative(U1,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U3 = U2 + dt/6 * computeFluxDerivative(U2,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U4 = U3 + dt/6 * computeFluxDerivative(U3,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U5 = U4 + dt/6 * computeFluxDerivative(U4,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
Ustar = 1/25 * U0 + 9/25 * U5;
U5 = 15 * Ustar - 5 * U5;
U6 = U5 + dt/6 * computeFluxDerivative(U5,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U7 = U6 + dt/6 * computeFluxDerivative(U6,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U8 = U7 + dt/6 * computeFluxDerivative(U7,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U9 = U8 + dt/6 * computeFluxDerivative(U8,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U = 6/10 * U9 + dt/10 * computeFluxDerivative(U9,flux,lam,gamma,dx,kappa,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn) + Ustar;
end
