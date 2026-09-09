function name = getSchemeName(flux_limiter,kappa)
% Decide scheme label based on flux limiter and kappa
% FIXED: body previously assigned to `Name` (capitalized), a different
% variable from the declared output `name`, so this function never
% returned a value.
if isscalar(flux_limiter)
    if flux_limiter ~= 18
        name = ['Roe-MUSCL_kappa=' num2str(kappa)];
    else
        name = 'Roe';
    end
elseif ~isscalar(flux_limiter) && any(flux_limiter==18)
    name = ['Roe_v_Roe-MUSCL_kappa=' num2str(kappa)];
else
    name = ['Roe-MUSCL_kappa=' num2str(kappa)];
end
end
