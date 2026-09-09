function U = SSPRKN(order,U,dt)
% NOTE: dead code, unreachable from runSolver.m (not in its time_int
% dispatch). Also broken as written: computeFluxDerivative is called
% with 1 argument instead of the 11 its real signature requires, the
% Uvec row-stacking doesn't match the 3xN conservative-variable
% convention used elsewhere, and the declared output `U` is never
% assigned (only local `Um` is computed). Left as-is / not ported.
alpha = SSPRKcoeffs(order);
Uvec = [U];
for i = 1:order
    Ui = Uvec(i,:) + dt/2 * computeFluxDerivative(Uvec(i,:));
    Uvec = [Uvec;Ui];
end
Um = zeros(size(U));
for i = 1:length(alpha)-1
    Um = Um + alpha(i)*Uvec(i,:);
end

Um = Um + alpha(end-1)*U;

end