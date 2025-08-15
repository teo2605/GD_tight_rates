function Ti = get_Ti(gL,gmu,i)
% gL = gamma * L; gmu = gamma * mu;
    if gmu == 0
        Ti = 2*i - (-1 + (1-gL).^(-2*i))/gL;
    else
        Ti = (-1 + (1-gmu).^(-2*i))/gmu - (-1 + (1-gL).^(-2*i))/gL ;
    end
end