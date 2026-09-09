function [URip12, ULip12] = weno3Rec(U)
% WENO3 interpolation scheme for linear advection equation du/dt+a.*du/dx=0
% Inputs:
%   U - advected quantity
%   i - cell face
% Outputs:
%   URip12 - reconstructed right cell face
%   ULip12 - reconstructed left cell face

%{
i-2 = 1:end-4
i-1 = 2:end-3
i = 3:end-2
i+1 = 4:end-1
i+2 = 5:end

URi = U(:, 3:end);      | Uip1, For the Roe Scheme
ULi = U(:, 2:end-1);    | Ui
URim1 = U(:, 2:end-1);  | Ui
ULim1 = U(:, 1:end-2);  | Uim1
%}

% Extract stencil points for all 3 variables at once
Uim2 = U(:,1:end-5);
Uim1 = U(:,2:end-4);
Ui = U(:,3:end-3);
Uip1 = U(:,4:end-2);
Uip2 = U(:,5:end-1);
Uip3 = U(:,6:end);

C1 = 1/3;
C2 = 2/3;
ep = 1e-6;
p = 2;

% FIXED: smoothness indicators and weights were previously computed only
% from component 1 (density) via Uip2(1,:) etc., then the SAME weights
% were applied to all 3 conservative variables - see weno5Rec.m for the
% full explanation (this is the same regression from the original
% WENO3_Rec.m, which computes smoothness/weights independently per
% conservative variable). Removing the "(1,:)" restriction below makes
% every formula operate on the full 3xN arrays (each row/variable
% independently), which is exactly per-component treatment and matches
% the original's per-component-unrolled result to floating-point
% precision.

% ===== RIGHT RECONSTRUCTION (at interface i+1/2) =====
IS1_R = (Uip2-Uip1).^2;
IS2_R = (Uip1-Ui).^2;

alpha1_R = C1./(ep+IS1_R).^p;
alpha2_R = C2./(ep+IS2_R).^p;
alphatot_R = alpha1_R + alpha2_R;

omega1_R = alpha1_R./alphatot_R;
omega2_R = alpha2_R./alphatot_R;

q1 = (3/2) .* Uip1 - (1/2) .* Uip2;
q2 = (1/2) .* Ui + (1/2) .* Uip1;
URip12 = omega1_R .* q1 + omega2_R .* q2;

% ===== LEFT RECONSTRUCTION (at interface i+1/2) =====
IS1_L = (Ui-Uim1).^2;
IS2_L = (Uip1-Ui).^2;

alpha1_L = C1./(ep+IS1_L).^p;
alpha2_L = C2./(ep+IS2_L).^p;
alphatot_L = alpha1_L + alpha2_L;

omega1_L = alpha1_L./alphatot_L;
omega2_L = alpha2_L./alphatot_L;

q1 = (3/2) .* Ui - (1/2) .* Uim1;
q2 = (1/2) .* Ui + (1/2) .* Uip1;
ULip12 = omega1_L .* q1 + omega2_L .* q2;

end
