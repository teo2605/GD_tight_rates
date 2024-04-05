function Ti = get_Ti(gL,gmu,i)
% k = mu/L; h = gamma*L;
    if gmu == 0
        Ti = (1 - (1-gL).^(-2*i))/gL + 2*i;
    else
        Ti = (1 - (1-gL).^(-2*i))/gL - (1 - (1-gmu).^(-2*i))/gmu;
    end
end