function fp = funPp(Pstar,P,rho,g)
% funPp - Derivative of pressure function for Riemann solver root finding.
%
% FIXED: the two branches below were swapped relative to funP.m's branch
% condition (funP uses `Pstar > P` for its shock branch; this used
% `P > Pstar`), and the shock-branch formula evaluated sqrt(Al/(Bl+P))
% and the denominator 2*(Bl+P) at P instead of at Pstar (Toro eq. 4.37
% requires evaluation at Pstar). Verified against a central-difference of
% funP.m.
%
% VECTORIZED: works elementwise on arrays (Pstar,P,rho all the same
% size, or scalars) as well as plain scalars, matching funP.m.

Al = 2./((g+1).*rho);
Bl = P.*(g-1)./(g+1);
al = sqrt(g.*P./rho);

fp = zeros(size(Pstar));
shock = Pstar > P;
fp(shock) = sqrt(Al(shock)./(Bl(shock) + Pstar(shock))) .* (1 - ((Pstar(shock) - P(shock))./(2*(Bl(shock) + Pstar(shock)))));
fp(~shock) = (1./(rho(~shock).*al(~shock))).*(Pstar(~shock)./P(~shock)).^(-(g+1)/(2*g));
end