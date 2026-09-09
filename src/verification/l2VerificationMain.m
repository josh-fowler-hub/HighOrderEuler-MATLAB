close all;
addpath(genpath(fullfile(pwd, 'src')));

cells = 1000;
L = 1;
gamma = 1.4;
ghost = 3;
cfl = 0.9;
flux = 1;
lam = 1;
kap = 1;
rec = 'WENO5';
CD_Term_Order = 1;
riemann_solver = 'Roe';

% FIXED: previously built xs as a bare interior-only linspace(0,1,cells),
% inconsistent with the ghost-padded xs convention used everywhere else
% in this codebase (setUpGrid.m). That now matters because the padded xs
% is needed to evaluate the MMS source term at the correct grid
% locations (see computeFluxDerivative.m's sourceTermFn).
[xsPadded, dx, n2, xidx] = setUpGrid(cells, ghost, 0, L);
xs = xsPadded(xidx); % interior-only x locations, for manufacturedSolutionSetUp
ys = zeros(size(xs));
zs = zeros(size(xs));

[W, W0] = manufacturedSolutionSetUp(L, xs, ys, zs);

Uexact = primToCons(gamma, W(1,:), W(2,:), W(5,:));
U = primToCons(gamma, W0(1,:), W0(2,:), W0(5,:));
[ULBC, URBC] = manufacturedSolutionBCs(U);
U = [ULBC, U, URBC];

% Compute exact solution with BCs once before the loop
[ULBC_exact, URBC_exact] = manufacturedSolutionBCs(Uexact);
Uexact_BCs = [ULBC_exact, Uexact, URBC_exact];

% FIXED: previously the manufactured-solution source term (sourceTerms.m)
% was never applied during time-stepping, so the numerical solution had
% no reason to track the manufactured "exact" solution at all - this
% wires it into computeFluxDerivative via the optional sourceTermFn
% argument now threaded through every time integrator.
sourceTermFn = @(x) sourceTerms(x);

steps = 0;
time = 0;
t = 0.0052;

while time < t
    steps = steps + 1;
    Uold = U;

    % Calculate the time step
    dt = calculateTimeStep(cfl, dx, time, t, U, gamma);
    time = time + dt;

    % Select the time integration method
    U = SSPRK65(U,flux,lam,gamma,dx,dt,kap,rec,CD_Term_Order,riemann_solver,xsPadded,time,sourceTermFn);

    % FIXED: previously concatenated fresh ghost columns onto U each
    % iteration (U = [ULBC, U, URBC]) without removing the existing
    % ones, so U grew by 6 columns every step instead of having its 3+3
    % ghost columns replaced in place (as applyBoundaryConditions.m does
    % for the main solver pipeline).
    [ULBC, URBC] = manufacturedSolutionBCs(U);
    U(:,1:ghost) = ULBC;
    U(:,end-ghost+1:end) = URBC;

    error = max(max(abs(U - Uold)));
    if error < 1e-3
        break
    end
end

plotVars(U,Uexact,xs);
