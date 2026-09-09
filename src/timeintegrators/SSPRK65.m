function U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn)
% SSPRK(6,5) strong stability preserving Runge-Kutta
% sourceTermFn - OPTIONAL, see computeFluxDerivative.m

if nargin < 13
    sourceTermFn = [];
end

% Maintains a global counter used for debugging/analysis
global time_int_counter;
time_int_counter = time_int_counter + 1;

U0 = U;
U1 = U0 + dt/2 *computeFluxDerivative(U0, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U2 = U1 + dt/2 *computeFluxDerivative(U1, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U3 = U2 + dt/2 *computeFluxDerivative(U2, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U4 = U3 + dt/2 *computeFluxDerivative(U3, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U5 = U4 + dt/2 *computeFluxDerivative(U4, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U6 = (1/9)*U0 + (2/5)*U1 + (4/9)*U3 +(2/45)*U5 + dt/45 *computeFluxDerivative(U5, flux, lam, gamma,dx,kap,rec,CD_Term_Order,riemann_solver,xs,time,sourceTermFn);
U = U6;
end
