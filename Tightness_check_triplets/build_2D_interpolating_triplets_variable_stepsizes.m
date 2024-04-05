function build_2D_interpolating_triplets_variable_stepsizes(L, mu, gamma, Delta, N)
% 2D interpolating triplets for variables stepsizes in [1, \bar{\gamma L}^1]
%% Default run
if nargin == 0
    L       =  3;         % upper curvature
    mu      = -1;         % lower curvature
    gammaL  = [1.1, 1.2, 1.3, 1.4, 1.7];        % variable (but fixed) stepsizes
    gamma   = gammaL / L;
    Delta   = 2;          % initial condition gap
    N       = 5;          % number of iterations
end
assert ( sum (gamma*L < 1) == 0 ); % step-sizes >= 1/L
%%
[~, wc_N] = get_wc_rate(L,mu,gamma,N,Delta);
U = sqrt(2*L*wc_N); % gradient norm
c_vals = ( 1 + (1-gamma*mu).*(1-gamma*L) ) ./ (2-gamma*L-gamma*mu);
theta_vals = acos(c_vals);
theta_ = @(i) theta_vals(i+1);
R = @(theta) ( [cos(theta) -sin(theta); ...
                sin(theta)  cos(theta)] ); % rotation matrix between gradients
%% build gradients
g_vals = zeros(2, N+1);
g_vals(:,1) = U * [1;0];
for i = 1 : N
    g_vals(:,i+1) = R ( (-1)^(i) * theta_(i-1) ) * g_vals(:,i);
end
g_ = @(i) (g_vals(:,i+1)); 
%% compute function values
p = gamma*L .* ( 2 + gamma*L .* gamma*mu ./ ( 1-gamma*mu - abs(1-gamma*L) ) );
p_ = @(i) p(i+1);
fN = 0;
fi_ = @(i) ( fN + wc_N * sum (p_(i:N-1)) );
f_vals = zeros(1,N+1);
for ii = 0 : N; f_vals(ii+1) = fi_(ii); end
f_ = @(i) (f_vals(:,i+1));
%% compute iterations x
x_vals = zeros(2, N+1);
for i = 1 : N;  x_vals(:,i+1) = x_vals(:,i) - gamma(i)*g_vals(:,i); end
x_vals = x_vals - x_vals(:,end); % make xN = 0
x_ = @(i) (x_vals(:,i+1));
%% triplets:
fprintf('Iterations x_i: \n'); disp(x_vals);
fprintf('Gradients g_i: \n'); disp(g_vals);
fprintf('Function values f_i: \n'); disp(f_vals);
%% Check interpolation conditions
L_ij = zeros(N+1);
R_ij = zeros(N+1);
% build interpolation inequalities, divided in lhs (L) and rhs (R)
for i = 0 : N
    for j = 0 : N
        L_ij(i+1, j+1) = f_(i)-f_(j) - g_(j).'* (x_(i) - x_(j));
        R_ij(i+1, j+1) = 1/(2*(1-mu/L)) * ...
                            ( 1/L    * norm(g_(i)-g_(j))^2 + ...
                              mu     * norm(x_(i)-x_(j))^2 - ...
                              2*mu/L * (g_(j)-g_(i)).'*(x_(j)-x_(i)) );
    end
end
Interp_conds = L_ij - R_ij;
check_interp_matrix(Interp_conds)
end
%%
function check_interp_matrix(I_check, tol)
    if nargin == 1;  tol = 1e-10;  end    
    num_violdated_interp_constr = sum(I_check < -tol, 'all');
    if num_violdated_interp_constr > 0
        fprintf("Violated interpolation conditions! = %d \n\n", num_violdated_interp_constr);
    else
        fprintf("\n \t\t All interpolation conditions verified! \n\n");
    end
end