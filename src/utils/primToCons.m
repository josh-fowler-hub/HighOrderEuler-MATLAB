function U = primToCons(gamma, rhos, us, Ps)
% FIXED: previously divided by rhos unconditionally, producing NaN at
% rho==0 (e.g. the vacuum region of an exact Riemann solution) even
% though the physically correct conservative state there is all-zero
% (no mass -> no momentum or energy).
es = zeros(size(rhos));
valid = rhos > 0;
es(valid) = Ps(valid)./(rhos(valid).*(gamma - 1));
U = [rhos; rhos.*us; rhos.*(0.5.*us.^2 + es)];
end
