function hs = compute_optimal_stepsize_hypoconvex(N, kappa, gL0)
% gL0 = initial guess
    if nargin == 0
        N = 100; 
        kappa = -0.10;
    end
    if nargin < 3
        h_start = roots([-kappa*(1+kappa), 3*kappa+(1+kappa)^2, -4*(1+kappa), 4]);
        h_start = h_start ( (h_start >=1) & (h_start <= 2/(1+max(0,kappa))) );
        h_start = unique(h_start); % to cover multiple equal solutions (convex case)
    else
        h_start = gL0;
    end
    L = 1;  mu = kappa*L;    
    f_to_min = @(h) ( -h*get_PN(L,mu,h/L,N) );
    
    opts.MaxIterations = 3000;
    opts.MaxFunctionEvaluations = 3000;
    opts.OptimalityTolerance = 1e-10;
    opts.Display = 'off';  
    
    hs = fmincon(f_to_min, h_start, [], [], [], [], 1, 2, [], opts);
end