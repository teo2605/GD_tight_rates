function PN = get_PN(L,mu,gamma,N)
if nargin < 4; N = 10; end        
if nargin < 3; gamma = 1.9; end
if nargin < 2; mu = -1; end
if nargin < 1; L = 1; end
%%
p = 0; PN = inf;
if gamma ~= 0
    p = ( 2 + gamma*L .* gamma*mu ./ ( 1-gamma*mu - abs(1-gamma*L) ) );
end
%% Rates for variable (but fixed) gamma
if length(unique(gamma)) > 1
    assert( length(gamma) == N );
    PN = sum ( (gamma*L) .* p );
else
    %% Rate for constant stepsizes gamma
    if mu <= 0 % nonconvex        
        if gamma*L <= 1   % short steps         
            PN = p * N; % p = (2 + gamma*L * gamma*mu / (gamma*L - gamma*mu) );
        else % p = (2-gamma*L) * (2-gamma*mu) / (2 - gamma*L - gamma*mu);            
            Nbar = get_Nbar( gamma*L, gamma*mu );  
            if N > Nbar
                % alternative formulas:
                PN_1 = p*(N - Nbar) + ...
                       (-1 + (1-gamma*L)^(-2*Nbar))/(gamma*L) + (1-(1-gamma*L)^2) / ( (1-gamma*mu)^2-(1-gamma*L)^2 ) * get_Ti(gamma*L, gamma*mu, Nbar);
                PN_2 = p*(N - Nbar - 1) + ...
                       (-1 + (1-gamma*L)^(-2*(Nbar+1)))/(gamma*L) + (-1+(1-gamma*L)^(-2)) / ((1-gamma*L)^(-2)-(1-gamma*mu)^(-2)) * get_Ti(gamma*L, gamma*mu, Nbar+1);
                PN = p* ( N - (-gamma*mu) * gamma*L / (gamma*L - gamma*mu) * sum( get_Ti(gamma*L, gamma*mu, (1:Nbar)) ) );
            else
                PN = (-1+(1-gamma*L).^(-2*N)) / (gamma*L);
            end
        end
    else % mu > 0; strongly convex
        PN = (-1 + (1-gamma*mu)^(-2*N)) / (gamma * mu);
    end
end
end