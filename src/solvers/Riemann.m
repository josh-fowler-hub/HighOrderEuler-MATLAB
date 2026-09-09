function [Pstar,ustar,aL,aR,Uexact] = Riemann(rhoL,uL,PL,rhoR,uR,PR,t,middle,R,g,cells,extra_prim_var_plot)

% Exact Riemann Solver for Euler Equations
% Inputs:
%   rhoL - left state density
%   uL - left state velocity
%   PL - left state pressure
%   rhoR - right state density
%   uR - right state velocity
%   PR - right state pressure
%   t - run time
%   middle - location of initial solution
%   R - gas constant
%   g - ratio of specific heats
%   cells - number of grid cells
%   extra_prim_var_plot - plot extra primitive fields

% Setup grid
x = linspace(0,1,cells);

aL = sqrt(g*PL/rhoL);
aR = sqrt(g*PR/rhoR);

% Solve the Riemann problem using a pressure root-finding procedure
% (see Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics").

du = uR - uL;
guess_type = 0;
Iter_Limit = 10000;
TOL = 1e-6;

% FIXED: was `pGuess` (lowercase g), an undefined function - only
% `Pguess` (capital P) is defined, so this call would error under
% case-sensitive MATLAB/Octave.
Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
[Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);

% FIXED: the original final fallback (guess_type >= 4) plotted the
% residual function and blocked on input() for a manually-typed pressure
% guess, which cannot run headless. That branch is replaced below with
% the existing robust bisection fallback (already used elsewhere in this
% codebase) followed by one Newton polish from its result, mirroring
% findPstar.m's escalation pattern. The loop is bounded by guess_type<4
% so it always terminates.
while (abs(err) > TOL || err < 0) && guess_type < 4
    if guess_type <= 2
        guess_type = guess_type + 1;
        Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
        [Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    else
        guess_type = guess_type + 1;
        [Pstar, err] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, 10*Iter_Limit, TOL);
        [Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    end
end

ustar = 0.5*(uL+uR)+0.5*(funP(Pstar,PR,rhoR,g)-funP(Pstar,PL,rhoL,g));

% FIXED: previously this just broadcast the left state across the whole
% domain (a stub - "this block can be extended if needed"). Now samples
% the true exact-Riemann-problem solution (shock/rarefaction fans + star
% region + contact, and the vacuum case) at each grid point via
% sampleExact.m (Toro Ch. 4), and converts to conservative variables so
% Uexact is directly comparable to the numerical solver's U (previously
% Uexact was left in primitive [rho,u,P] form despite U being
% conservative [rho, rho*u, E] - the L2Norm/LinfNorm calls on these two
% arrays would otherwise compare mismatched quantities).
S = (x - middle) / t;
[rho, u, p] = sampleExact(S, rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pstar, ustar, g);
Uexact = primToCons(g, rho, u, p);

if extra_prim_var_plot
    % Additional plotting is handled outside this function
end
end
