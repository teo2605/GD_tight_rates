function wc = gradient_descent_dynamic_stepsizes_PEP(L, mu, gamma, Delta, N, use_opt)
%% Default run
if nargin == 0
    close all; clc;
    L       =  1;         % upper curvature
    mu      = -1e-3;      % lower curvature
    gamma   = 1;          % constant step-size (for consistency; not used like this for N>1)
    Delta   = 1;          % initial condition gap
    N       = 20;         % number of iterations
    use_opt = true;       % if false: init. cond. (f0-fN)<=Delta; if true: (f0-f*)<=Delta
end
if length(gamma) == 1 % if ct stepsize we switch to dynamic sequence
    kappa = mu/L;
    s_vec = get_sequence_stepsizes (kappa, N);
    gamma_seq = s_vec(2:end); % s_vec(1) = 0;
    
    if kappa < 0 % compute asymptotically optimal stepsize
        hs = roots([-kappa*(1+kappa), 3*kappa+(1+kappa)^2, -4*(1+kappa), 4]);
        hs = unique( hs ( (hs >=1) & (hs <= 2/(1+max(0,kappa))) ) );
    else
        hs = 2/(1+kappa); % the limit of the sequence for (strongly) convex funtions
    end    
    gamma = min(gamma_seq, hs);
end
gamma_vec = gamma;
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
P.TraceHeuristic(0); % 1 or 2 to get a low dimensional solution
result = P.solve(2);
wc = result.WCperformance;
%% (7) Evaluate the output
if strcmp(result.solverDetails.info, 'Unbounded objective function (MOSEK)')
    if ~isempty(gamma); gamma = gamma(1); end
    fprintf("\n\n \t\t UNBOUNDED -- Setup: gamma = %.4f, mu1 = %.4f \n\n", gamma, mu);
end
%% Compare with theoretical rate
[wc_opt, wc_N] = get_wc_rate(L,mu,gamma,N,Delta);
if use_opt
    wc_th = wc_opt;
else
    wc_th = wc_N;
end
fprintf( "\t\t PEP rate: \t\t\t  wc =  %.7f \n" , wc );
fprintf( "\t\t Theoretical rate:\t  wc_th =  %.7f \n" , wc_th );
fprintf( "\t\t Error: \t\t\t  wc-wc_th = %.7e \n" , wc - wc_th );

end