function fluxLim = getFluxName(flux)
% Map flux-limiter ID to a string name
% FIXED: body previously assigned to `flux_lim`, a different variable from
% the declared output `fluxLim`, so this function never returned a value.
switch flux
    case 1
        fluxLim = 'MM';
    case 2
        fluxLim = 'VL';
    case 3
        fluxLim = 'BJ';
    case 4
        fluxLim = 'SB';
    case 5
        fluxLim = 'VA2';
    case 6
        fluxLim = 'VA1';
    case 7
        fluxLim = 'CHARM';
    case 8
        fluxLim = 'HCUS';
    case 9
        fluxLim = 'HQUICK';
    case 10
        fluxLim = 'Koren';
    case 11
        fluxLim = 'MC';
    case 12
        fluxLim = 'Osher';
    case 13
        fluxLim = 'Ospre';
    case 14
        fluxLim = 'Smart';
    case 15
        fluxLim = 'Sweby';
    case 16
        fluxLim = 'UMIST';
    case 17
        fluxLim = 'GMM';
    case 18
        fluxLim = '';
    otherwise
        fluxLim = '';
end
end
