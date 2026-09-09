function [rho,u,p] = sampleExact(S, rhoL, uL, PL, aL, rhoR, uR, PR, aR, Pstar, ustar, g)
% sampleExact - Samples the exact solution of the 1D Riemann problem
% (Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics",
% Ch. 4) at similarity variable S = (x - x0)/t, given the already-solved
% star-region pressure/velocity (Pstar, ustar). S may be a vector; the
% outputs are the same shape.
%
% Handles all four wave-pattern combinations (shock or rarefaction on
% each side of the contact) as well as the vacuum-generation case
% (Toro Sec. 4.8, when the two states separate fast enough that a vacuum
% forms between them).

S = S(:).';
rho = zeros(size(S));
u = zeros(size(S));
p = zeros(size(S));

vacuum = (uR - uL) >= (2*aL/(g-1) + 2*aR/(g-1));

if vacuum
    % Toro Sec. 4.8: left/right rarefaction fans expand into a vacuum
    % region rather than meeting at a star region.
    SLstar = uL - 2*aL/(g-1);
    SRstar = uR + 2*aR/(g-1);

    mLdata = S < (uL - aL);
    mLfan  = (~mLdata) & (S < SLstar);
    mVac   = (S >= SLstar) & (S <= SRstar);
    mRfan  = (~mVac) & (S <= (uR + aR)) & (S > SRstar);
    mRdata = S > (uR + aR);

    rho(mLdata) = rhoL; u(mLdata) = uL; p(mLdata) = PL;

    c = (2/(g+1)).*(aL + 0.5*(g-1).*(uL - S(mLfan)));
    u(mLfan) = (2/(g+1)).*(aL + 0.5*(g-1)*uL + S(mLfan));
    rho(mLfan) = rhoL.*(c./aL).^(2/(g-1));
    p(mLfan) = PL.*(c./aL).^(2*g/(g-1));

    rho(mVac) = 0; p(mVac) = 0; u(mVac) = 0;

    c = (2/(g+1)).*(aR - 0.5*(g-1).*(uR - S(mRfan)));
    u(mRfan) = (2/(g+1)).*(-aR + 0.5*(g-1)*uR + S(mRfan));
    rho(mRfan) = rhoR.*(c./aR).^(2/(g-1));
    p(mRfan) = PR.*(c./aR).^(2*g/(g-1));

    rho(mRdata) = rhoR; u(mRdata) = uR; p(mRdata) = PR;
    return
end

leftShock = Pstar > PL;
rightShock = Pstar > PR;

% ---- Left of the contact (S <= ustar) ----
mLeft = S <= ustar;

if leftShock
    SL = uL - aL*sqrt((g+1)/(2*g)*(Pstar/PL) + (g-1)/(2*g));
    mLdata = mLeft & (S <= SL);
    mLstar = mLeft & ~(S <= SL);

    rho(mLdata) = rhoL; u(mLdata) = uL; p(mLdata) = PL;

    ratio = Pstar/PL;
    rhoLstar = rhoL*(ratio + (g-1)/(g+1)) / ((g-1)/(g+1)*ratio + 1);
    rho(mLstar) = rhoLstar; u(mLstar) = ustar; p(mLstar) = Pstar;
else
    SHL = uL - aL;
    astarL = aL*(Pstar/PL)^((g-1)/(2*g));
    STL = ustar - astarL;

    mLdata = mLeft & (S <= SHL);
    mLfan  = mLeft & (S > SHL) & (S < STL);
    mLstar = mLeft & (S >= STL);

    rho(mLdata) = rhoL; u(mLdata) = uL; p(mLdata) = PL;

    c = (2/(g+1)).*(aL + 0.5*(g-1).*(uL - S(mLfan)));
    u(mLfan) = (2/(g+1)).*(aL + 0.5*(g-1)*uL + S(mLfan));
    rho(mLfan) = rhoL.*(c./aL).^(2/(g-1));
    p(mLfan) = PL.*(c./aL).^(2*g/(g-1));

    rhoLstar = rhoL*(Pstar/PL)^(1/g);
    rho(mLstar) = rhoLstar; u(mLstar) = ustar; p(mLstar) = Pstar;
end

% ---- Right of the contact (S > ustar) ----
mRight = S > ustar;

if rightShock
    SR = uR + aR*sqrt((g+1)/(2*g)*(Pstar/PR) + (g-1)/(2*g));
    mRdata = mRight & (S >= SR);
    mRstar = mRight & ~(S >= SR);

    rho(mRdata) = rhoR; u(mRdata) = uR; p(mRdata) = PR;

    ratio = Pstar/PR;
    rhoRstar = rhoR*(ratio + (g-1)/(g+1)) / ((g-1)/(g+1)*ratio + 1);
    rho(mRstar) = rhoRstar; u(mRstar) = ustar; p(mRstar) = Pstar;
else
    SHR = uR + aR;
    astarR = aR*(Pstar/PR)^((g-1)/(2*g));
    STR = ustar + astarR;

    mRdata = mRight & (S >= SHR);
    mRfan  = mRight & (S < SHR) & (S > STR);
    mRstar = mRight & (S <= STR);

    rho(mRdata) = rhoR; u(mRdata) = uR; p(mRdata) = PR;

    c = (2/(g+1)).*(aR - 0.5*(g-1).*(uR - S(mRfan)));
    u(mRfan) = (2/(g+1)).*(-aR + 0.5*(g-1)*uR + S(mRfan));
    rho(mRfan) = rhoR.*(c./aR).^(2/(g-1));
    p(mRfan) = PR.*(c./aR).^(2*g/(g-1));

    rhoRstar = rhoR*(Pstar/PR)^(1/g);
    rho(mRstar) = rhoRstar; u(mRstar) = ustar; p(mRstar) = Pstar;
end
end
