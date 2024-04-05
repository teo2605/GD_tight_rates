function wc = gradient_descent_PEP(L, mu, gamma, Delta, N, use_opt)
% Run PEP over gradient descent for any curvatures L and mu and any
% stepsize. It returns:
% $$ wc = \min_{0\leq i \leq N} \frac{\{\nabla f(x_i)\}^2}{2L}$$
% and compares it with the theoretical rate obtained for the cases:
%   1. constant stepsizes $\gamma$ with $\gamma L \in (0,2)$
%   2. variable stepsizes $\gamma_i$ with $\gamma_i L \in (0,\bar{\gamma L}_i]$, when $\mu \leq 0$.
% Input parameter use_opt distinguishes between: 
%   (i) true: using optimal $f_*$ in computations (tightly);
%   (ii) false: using the gap $f(x_0)-f(x_N)$
%
    %% Default run
    if nargin == 0
        close all; clc;
        L       =  1;         % upper curvature
        mu      = -1;         % lower curvature
        gamma   = 1.85;       % step-size
        Delta   = 1;          % initial condition gap
        N       = 10;         % number of iterations
        use_opt = true;       % if false: init. cond. (f0-fN)<=Delta; if true: (f0-f*)<=Delta
    end
    %%
    if length(gamma) == 1 
        gamma_vec = gamma * ones(1, N);
    else
        gamma_vec = gamma; % to use for variable step-sizes
    end
    %% (0) Initialize an empty PEP
    P = pep();
    % (1) Set up the objective function
    param.L  = L;     % Smoothness parameter
    param.mu = mu;    % Curvature parameter
    % f is the objective function
    f = P.DeclareFunction('Hypoconvex', param); % alternatively, strongly-convex but without condition mu>0
    % (2) Declare data structures
    x_       = cell(N+1,1);
    g_       = cell(N+1,1);
    f_       = cell(N+1,1);
    if use_opt
        % (3.1) Define an optimal point
        [~, f_s]= f.OptimalPoint();
    end
    % (3.2) Define the starting point
    x_str = @(i)(sprintf('x_%d', i-1));
    x_{1}         = P.StartingPoint();		 % x0=x_{1} is some starting point
    [g_{1},f_{1}] = f.oracle(x_{1}, x_str(1));
    if use_opt;  P.AddConstraint( f_{1} - f_s - 1/(2*L) * g_{1}^2 >= 0 ); end % optimality constraint
    % (4) Algorithm -- Iterations of GM
    for i = 1 : N
        x_{i+1}            = x_{i} - gamma_vec(i) * g_{i}; % GM update
        [g_{i+1}, f_{i+1}] = f.oracle(x_{i+1}, x_str(i+1));    
        if use_opt; P.AddConstraint( f_{i+1} - f_s - 1/(2*L) * g_{i+1}^2 >= 0 ); end % optimality constraint    
    end
    % (5.1) Set up the performance meassure
    for i = 1 : N+1
        P.PerformanceMetric( g_{i}^2 / (2*L) );
    end
    % (5.2) Set up the initial condition
    if use_opt
        P.InitialCondition( f_{1} - f_s <= Delta );
    else
        P.InitialCondition( f_{1} - f_{N+1} <= Delta );
    end
    %% (6) Solve the PEP
    P.TraceHeuristic(0); % force low-dimensional solution (if 1 or 2); 
    result = P.solve(2);
    wc = result.WCperformance;
    %% (7) Evaluate the output
    if strcmp(result.solverDetails.info, 'Unbounded objective function (MOSEK)') % if using MOSEK as SDP solver
        if ~isempty(gamma); gamma = gamma(1); end
        fprintf("\n\n \t\t UNBOUNDED -- Setup: gamma = %f, mu1 = %f \n\n", gamma, mu);
    end
    %% Compare with theoretical rate
    [wc_th_opt, wc_th_N] = get_wc_rate(L,mu,gamma,N,Delta);
    if use_opt
        wc_th = wc_th_opt;
    else
        wc_th = wc_th_N;
    end
    fprintf( "\t\t PEP rate: \t\t\t  wc =  %.7f \n" , wc );
    fprintf( "\t\t Theoretical rate:\t  wc_th =  %.7f \n" , wc_th );
    fprintf( "\t\t Error: \t\t\t  wc-wc_th =  %.7e \n" , wc - wc_th );
end