clear; clc; close all
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% setup
kappa = -.001; L = 1; Delta = 1; mu = kappa*L;
Nvec = [1,2,5,10,20,30,40,50,100];
% Nvec = 1 : 70;
Nmax = max(Nvec);
%% weights of gradient norms after one iteration
sigma_0 = @(h) ( h.*( (2-h).*(2-kappa*h) - 1) ./ (2-h*(1+kappa)) );
sigma_1 = @(h) (h./(2-h*(1+kappa)) );
p = @(h) (sigma_0(h) + sigma_1(h));
%% get stepsizes
h_var_vec = zeros(1,Nmax+1);
for N = 1 : Nmax
    h_var_vec(N+1) = get_next_stepsize(h_var_vec(N), kappa);
end
%% compute the asymptotically optimal stepsize
hs = roots([-kappa*(1+kappa), 3*kappa+(1+kappa)^2, -4*(1+kappa), 4]);
hs = hs ( (hs >=1) & (hs <= 2/(1+max(0,kappa))) );
hs = unique(hs); % to cover multiple equal solutions (convex case)

s = @(i) ( h_var_vec(i+2) ); % dynamic stepsize sequence
Ntilde = find(h_var_vec > hs, 1) - 3; % iteration with the switch in stepsize policy
Nvec = sort(unique([Ntilde+1, Ntilde+2, Nvec]));
Nvec_sel = Nvec + 1;
Nmax = max(Nvec);
%% worst-cases
Nvec = 0 : Nmax; % for plots we compute for all iterations
P_den_hs = 0*Nvec;
P_den_h_dyn = 0*Nvec;
P_den_seq = 0*Nvec;
hs_N = 0 * Nvec;
P_den_hs_N = 0 * Nvec;
for NN = 1 : length(Nvec)
    % clc
    N = Nvec(NN);
    P_den_hs(NN)  = 1 + hs*get_PN(L,mu,hs/L,N); % 
    P_den_seq(NN) = 1 + sigma_1(h_var_vec(NN));
    P_den_h_dyn(NN) = 1 + max ( sigma_1(h_var_vec(NN)), p(hs)*(N-Ntilde-1) + sigma_1(s(Ntilde)) );    
    % opt ct stepsize: hs = hs(N,kappa)
    hs_N(NN) = compute_optimal_stepsize_hypoconvex(N,kappa);
    P_den_hs_N(NN) = 1 + hs_N(NN)*get_PN(L,mu,hs_N(NN)/L,N);
end
%% Compare stepsizes
figure(1); clf;  leg_names = {};
h_var_vec = min(2, h_var_vec);
plot( (0:Nmax), hs*ones(1,Nmax+1), '-bo', 'LineWidth', 1); hold on; leg_names{end+1} = 'asympt. $(\gamma L)_*(\kappa)$';
plot( (0:Nmax), hs_N, '-mo', 'LineWidth', 1); leg_names{end+1} = 'ct. $(\gamma L)_*(N,\kappa)$';
plot( (0:Nmax), h_var_vec, '-k*', 'LineWidth', 1); leg_names{end+1} = '$s_{N-1}(\kappa)$';
plot( (0:Nmax), min(h_var_vec, hs*ones(1,Nmax+1)), '-r*', 'LineWidth', 1); hold on; leg_names{end+1} = 'dynamic';
legend(leg_names, 'Location','southeast')
xlabel('$N$', 'Interpreter','latex')
ylabel('Stepsize $(\gamma L)_N$', 'Interpreter','latex')
title_str = sprintf("Sequence of stepsizes for $\\kappa=%.e$", kappa);
title(title_str, 'Interpreter','latex')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 2;
grid on
xlim([1 Nmax]); %ylim([1.5 2])
xticks_init = xticks; 
xticks([1, xticks_init])
%% Compare denomniators (semilogy)
figure(2); clf; leg_names = {};
semilogy(Nvec, 1./P_den_hs, 'b', 'LineWidth', 2); hold on; leg_names{end+1} = 'asympt. $(\gamma L)_*(\kappa)$';
semilogy(Nvec, 1./P_den_hs_N,'-m', 'LineWidth',2); leg_names{end+1} = "ct. $(\gamma L)_*(N,\kappa)$";
semilogy(Nvec, 1./P_den_seq,'k', 'LineWidth',2); leg_names{end+1} = "$s_{N-1}(\kappa)$";
semilogy(Nvec, 1./P_den_h_dyn, 'r--', 'LineWidth', 2); hold on; leg_names{end+1} = 'dynamic';
legend(leg_names, 'Interpreter','latex','FontSize',12, 'Location','northeast')
xlabel('Number of iterations $N$')
ylabel('$(\gamma L) \, P_N(\gamma L, \kappa)$')
title_str = sprintf("Denominators of rates ($\\kappa=%.e$)", kappa);
title(title_str, 'Interpreter','latex')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 2;
grid on
%% Table: compare denominators
Iterations = num2str(Nvec_sel' - 1);
var_stepsizes = ( min( h_var_vec(Nvec_sel), hs ) )';
Den_var_steps = P_den_h_dyn(Nvec_sel)';
wc_var_stepsizes_N_its = 1 ./ Den_var_steps;

ct_N_stepsizes = hs_N(Nvec_sel)';
Den_ct_N_steps = P_den_hs_N(Nvec_sel)';

ct_stepsizes = hs + 0*(Nvec_sel)';
Den_ct_steps = P_den_hs(Nvec_sel)';

Ratio = Den_var_steps ./ Den_ct_N_steps * 100;
rate_comparison = table( Iterations, ...
    ct_stepsizes, Den_ct_steps, ...
    ct_N_stepsizes, Den_ct_N_steps, ...
    var_stepsizes, Den_var_steps, ...
    Ratio, ...
    wc_var_stepsizes_N_its);
n_decimal = 3;
table_compare_rates = varfun(@(x) num2str(x, ['%' sprintf('.%df', n_decimal)]), rate_comparison);
table_compare_rates.Properties.VariableNames = rate_comparison.Properties.VariableNames;
table_compare_rates.Properties.RowNames = rate_comparison.Properties.RowNames;
disp(table_compare_rates)