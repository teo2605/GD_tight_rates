function build_3D_interpolating_triplets(L, mu, gamma, Delta, N)
%% Default run
if nargin == 0
    clc;
    L       =  2;         % upper curvature
    mu      = -.5;         % lower curvature
    gammaL  = 1.93;        % (normalized) step-size
    gamma   = gammaL / L; % step-size
    Delta   = 3;          % initial condition gap
    N       = get_Nbar(gamma*L, gamma*mu) + 5;         % number of iterations
end
%%
Nbar = get_Nbar(gamma*L, gamma*mu);
n = 1-gamma*L; m = 1-gamma*mu;
[~, wc] = get_wc_rate(L,mu,gamma,N,Delta);
U = sqrt(2*L*wc);
Utilde = sqrt( (m^2-1) * (1-n^2) / (m^2-n^2) );
%% define first Nbar+1 gradients
g_quad_ = @(i) U* Utilde * (  [ n.^(i-Nbar) / sqrt(1-n^2) ; ...
                                m.^(i-Nbar) / sqrt(m^2-1) ; ...
                                0] );
g_vals = zeros(3,N+1);
for i = 0 : Nbar+1;  g_vals(:,i+1) = g_quad_(i); end
%% define the last (N-Nbar) gradients
q_ = @(i) (-1)^(i-(Nbar+1)) * sqrt( (m^2-1)*(1-n^(2*(i-(Nbar+1))))/(1-n^2) ) ;
R_ = @(i)  1/( 1+ q_(i)^2 ) * blkdiag(0, ...
                                        [ 1,   -q_(i) ; ...
                                          q_(i) ,   1] ) ;
for i = Nbar+1 : N
    g_vals(:,i+1) = ( n*eye(3) + (m-n)*R_(i) ) * g_vals(:,i);
end
g_ = @(i) (g_vals(:,i+1)); 
%% function values
Ti = @(i) ( get_Ti(gamma*L, gamma*mu, i) );
p = gamma*L * (2-gamma*L) * (2-gamma*mu) / (2-gamma*L-gamma*mu);
fN = 0;
fi_ = @(i) ( fN + p*U^2/(2*L) * ( (N-i) + (gamma*L * gamma*mu)/(gamma*L-gamma*mu) * sum ( Ti(1:Nbar-i)  ) ) );
f_vals = zeros(1,N+1);
for ii = 0 : N;  f_vals(ii+1) = fi_(ii); end
f_ = @(i) (f_vals(:,i+1)); 
%% iterations
x_vals = zeros(3, N+1);
for i = 1 : N;  x_vals(:,i+1) = x_vals(:,i) - gamma*g_vals(:,i); end
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
    if num_violdated_interp_constr > 0 || sum(isnan(I_check), 'all')
        fprintf("Violated interpolation conditions! = %d \n\n", num_violdated_interp_constr);
    else
        fprintf("\n \t\t All interpolation conditions verified! \n\n");
    end
end