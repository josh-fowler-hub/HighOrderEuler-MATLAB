function [Pstar_root,error,Iter] = newtonRaphsonP(Pstar0,PL,PR,rhoL,rhoR,g,du,iter_limit,tol)
% newtonRaphsonP - Newton-Raphson root-finder for Riemann problem pressure.
%
% VECTORIZED: Pstar0,PL,PR,rhoL,rhoR,du may be scalars or arrays of the
% same size. Each element is iterated independently, exactly as if
% newtonRaphsonP had been called once per element in a scalar loop (each
% element's update at iteration k depends only on its own previous
% value, not on other elements), but with each outer while-iteration
% updating every not-yet-converged element in one vectorized pass
% instead of one scalar Octave function call per element per iteration -
% this is what makes adaptiveRiemannSolver.m's per-interface Newton
% solve fast instead of the dominant cost of a full run.
% Previously-converged elements are frozen (not further updated), so
% per-element results are identical to the original scalar-loop version.

Pstar = Pstar0;
error = 2*tol*ones(size(Pstar));
converged = false(size(Pstar));
Iter = 0;

while Iter < iter_limit && ~all(converged)
    Iter = Iter + 1;
    active = ~converged;

    f = funP(Pstar(active),PL(active),rhoL(active),g)+funP(Pstar(active),PR(active),rhoR(active),g)+du(active);
    fp = funPp(Pstar(active),PL(active),rhoL(active),g)+funPp(Pstar(active),PR(active),rhoR(active),g);

    Pstarkm1 = Pstar(active);
    Pstarnew = Pstarkm1 - f./fp;
    error(active) = abs((Pstarnew-Pstarkm1)./(0.5*(Pstarnew+Pstarkm1)));
    Pstar(active) = Pstarnew;
    converged(active) = error(active) <= tol;
end
Pstar_root = Pstar;
end