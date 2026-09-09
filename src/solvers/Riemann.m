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

<<<<<<< HEAD
% FIXED: was `pGuess` (lowercase g), an undefined function - only
% `Pguess` (capital P) is defined, so this call would error under
% case-sensitive MATLAB/Octave.
Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
[Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
=======
Pstar = Pguess(guess_type, PL, PR, uL, uR, rhoL, rhoR, g, du);
[Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
>>>>>>> 0044426726c1b608e9272ae189e265c67a0d550d

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
<<<<<<< HEAD
        [Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    else
        guess_type = guess_type + 1;
        [Pstar, err] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, 10*Iter_Limit, TOL);
        [Pstar, err] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
=======
        [Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
    elseif guess_type > 2 && guess_type < 4
        guess_type = guess_type + 1;
        disp('Newton''s Method failed... Attempting Bisection Method.');
        [Pstar, error] = bisectionP(PL, PR, uL, uR, rhoL, rhoR, g, du, 10*Iter_Limit, TOL);
    else
        disp('Cannot Numerically find Pstar, Please Solve Graphically');
        Pstar0 = Pguess(0, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar1 = Pguess(1, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar2 = Pguess(2, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar3 = Pguess(3, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstar4 = Pguess(4, PL, PR, uL, uR, rhoL, rhoR, g, du);
        Pstara = min([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
        Pstarb = max([Pstar0,Pstar1,Pstar2,Pstar3,Pstar4]);
        Pstars = linspace(Pstara, Pstarb, 1000);
        ErrFunc = zeros(size(Pstars));
        for i = 1:length(Pstars)
            ErrFunc(i) = funP(Pstars(i),PL,rhoL,g)+funP(Pstars(i),PR,rhoR,g)+du;
        end
        figure(5)
        plot(Pstars, ErrFunc);
        grid off;
        drawnow;
        Pstar = input('Pstar0 Guess from Graph: ');
        [Pstar, error] = newtonRaphsonP(Pstar,PL,PR,rhoL,rhoR,g,du,Iter_Limit,TOL);
        break;
>>>>>>> 0044426726c1b608e9272ae189e265c67a0d550d
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
