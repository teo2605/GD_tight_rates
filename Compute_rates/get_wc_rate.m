function [wc_star, wc_N, DS, DN] = get_wc_rate(L,mu,gamma,N,Delta)
    % compute worst-case rate 
    %   $\min_i\{\|\nabla f(x_i)\|^2\}/(2L)$ 
    % under initial conditions:
    %   $f(x_0) - f_* \leq \Delta$ or 
    %   $f(x_0) - f(x_N) \leq \Delta$
    %%
    if nargin < 5; Delta = 1; end
    if nargin < 4; N = 5; end        
    if nargin < 3; gL = [1.1, 1.2, 1.3, 1.4, 1.7]; end
    if nargin < 2; mu = -1; end
    if nargin < 1; L = 1; end
    if nargin < 3; gamma = gL / L; end
    %%
    PN = get_PN(L,mu,gamma,N);
    if length(unique(gamma)) > 1 % denominator for variable stepsizes
        DN = PN;
    else % denominator for constant stepsizes
        DN = min ( gamma*L * PN, -1+(1-gamma*L).^(-2*N) ); % criterion f0-fN
    end
    DS = 1 + DN; % criterion f0-f*
    wc_N = Delta ./ DN; 
    wc_star = Delta ./ DS;    
end