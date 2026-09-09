function [Fip12, Fim12] = hllSolver(U,URi,ULi,URim1,ULim1,xs,g,t)
% hllSolver - HLL approximate Riemann solver (Toro Ch. 10), computed
% independently at both the i+1/2 interface (from ULi/URi) and the
% i-1/2 interface (from ULim1/URim1).
%
% FIXED: previously allocated Fip12/Fim12 as zeros(size(U)) - the FULL
% padded array width (including ghost cells) - then wrote into them
% using maskL/maskR/maskM, which are sized to match the RECONSTRUCTED
% (interior-only, narrower) URi/ULi arrays. Because the mask was
% narrower than U's width, the logical-indexed assignment silently wrote
% the correct flux values into the wrong (shifted) columns and left the
% last few columns permanently zero - corrupting the flux array with a
% column misalignment that caused near-instant numerical blow-up. Also
% previously used one shared SL/SR pair (from adaptiveRiemannSolver
% spanning Lim1..Ri) for both the i+1/2 and i-1/2 fluxes, and called
% adaptiveRiemannSolver expecting 4 outputs from its 3-output signature
% (same bug as noted in the original hllcSolver.m). Restructured to the
% same per-interface pattern now used by hllcSolver.m/osherSolver.m -
% Fip12/Fim12 sized to match URi/ULi directly, with wave speeds computed
% independently at each interface - which avoids the misalignment bug
% entirely and keeps the three approximate solvers consistent. The dead
% full-U stencil computation (Uim3..Uip3, Fim3..Fip3), never used in the
% original either, has been removed.
%
% Inputs:
%   U      - full conservative state array (unused, kept for signature
%            uniformity with the other Riemann solvers)
%   URi    - right reconstructed state at interface i+1/2
%   ULi    - left reconstructed state at interface i+1/2
%   URim1  - right reconstructed state at interface i-1/2
%   ULim1  - left reconstructed state at interface i-1/2
%   xs     - grid locations (unused)
%   g      - ratio of specific heats
%   t      - current time (unused)
%
% Outputs:
%   Fip12, Fim12 - numerical fluxes on the right and left side of each cell

Fip12 = hllFlux(ULi, URi, g);
Fim12 = hllFlux(ULim1, URim1, g);
end

function F = hllFlux(UL, UR, g)
[rhoL, uL, pL, aL, ~, ~, ~, ~] = consToPrim(UL, g);
[rhoR, uR, pR, aR, ~, ~, ~, ~] = consToPrim(UR, g);

% FIXED: previously used SL=ustar-astarL, SR=ustar+astarR directly, with
% no floor/ceiling against the raw data wave speeds uL-aL/uR+aR. Since
% the adaptive (PVRS/TRRS/TSRS) estimate is only approximate, it can
% under-bound the true wave fan at some interfaces (confirmed
% numerically: ~3% of interfaces on the Sod problem), which selects the
% wrong flux branch there and destabilizes the scheme within a handful
% of timesteps. This is the standard Toro/Einfeldt robust HLL wave-speed
% estimate (Toro Ch. 10) - the true outer bound is never weaker than the
% raw left/right data speeds.
[ustar, astarL, astarR] = adaptiveRiemannSolver(rhoL, rhoR, uL, uR, aL, aR, pL, pR, g);
SL = min(uL - aL, ustar - astarL);
SR = max(uR + aR, ustar + astarR);

FL = Fvec(rhoL, uL, pL, g);
FR = Fvec(rhoR, uR, pR, g);

maskL = SL >= 0;
maskR = SR <= 0;
maskM = ~(maskL | maskR);

F = zeros(size(UL));
F(:,maskL) = FL(:,maskL);
F(:,maskR) = FR(:,maskR);
F(:,maskM) = (SR(maskM).*FL(:,maskM) - SL(maskM).*FR(:,maskM) + SL(maskM).*SR(maskM).*(UR(:,maskM) - UL(:,maskM))) ./ (SR(maskM) - SL(maskM));
end
