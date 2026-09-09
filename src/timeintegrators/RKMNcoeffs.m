function alpha = RKMNcoeffs(m)
    % NOTE: dead code, not called anywhere in src/ or main.m. Also the
    % outer `for t = 1:1000` loop below is a no-op (its body doesn't
    % depend on t, so it recomputes the same result 1000 times). Left
    % as-is / not ported.
    % alphamk = 1/k * amm1km1;
    % alphammm1 = 1/factorial(m);
    % am0 = 1 - sum(amk);
    alpha = zeros([m,m]);
    alpha(1,1) = 0;
    alpha(1,2) = 1;
    for t = 1:1000
        for i = 2:m
            for j = 2:m-2
                alpha(i,j) = (1/j)*alpha(i-1,j-1);
                alpha(m,j) = (2/j)*alpha(m-1,j-1);
            end
        end
    end
    alpha(m,m-1) = (2/m)*alpha(m-1,m-2);
    for i = 1:m-1
        alpha(i,1) = 1 - sum(alpha(i,2:end));
    end
end