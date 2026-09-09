function f = funP(Pstar,P,rho,g)
% funP - Pressure function for Riemann solver root finding.
%
% VECTORIZED: works elementwise on arrays (Pstar,P,rho all the same
% size, or scalars) as well as plain scalars - the scalar case is just
% the N=1 case of the array case, so existing scalar call sites are
% unaffected. Enables vectorized Newton iteration (see
% newtonRaphsonP.m) instead of per-index scalar loops.

Ar = 2./((g+1).*rho);
Br = P.*(g-1)./(g+1);
ar = sqrt(g.*P./rho);

f = zeros(size(Pstar));
shock = Pstar > P;
f(shock) = (Pstar(shock)-P(shock)).*sqrt(Ar(shock)./(Pstar(shock)+Br(shock)));
f(~shock) = ((2.*ar(~shock))./(g-1)).*((Pstar(~shock)./P(~shock)).^((g-1)/(2*g)) - 1);
end