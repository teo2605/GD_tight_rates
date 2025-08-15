clear; clc; close all
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% setup
kappa = -1e-3; L = 1; Delta = 1; mu = kappa*L;
Nvec = [1,2,5,10,20,30,40,50,100,300];
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
plot( (0:Nmax), hs*ones(1,Nmax+1), 'm', 'LineWidth', 4); hold on; leg_names{end+1} = 'Asymptotically $(\gamma L)_*(\kappa)$';
plot( (0:Nmax), hs_N, '-b', 'LineWidth', 4); leg_names{end+1} = 'Best constant $(\gamma L)_*(N,\kappa)$';
plot( (0:Nmax), h_var_vec, 'k', 'LineWidth', 4); leg_names{end+1} = '$s_{N}(\kappa)$';
plot( (0:Nmax), min(h_var_vec, hs*ones(1,Nmax+1)), 'r:', 'LineWidth', 4); hold on; leg_names{end+1} = 'Scheduled';
legend(leg_names, 'Location','southeast', 'fontsize', 17)
xlabel('Number of iterations $N$', 'Interpreter','latex')
ylabel('Stepsize $\gamma L$', 'Interpreter','latex')
title_str = sprintf("Sequence of stepsizes for $\\kappa=%.e$", kappa);
% title(title_str, 'Interpreter','latex')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 17;
% ax.GridLineWidth = 1;
grid on
xlim([1 Nmax]); %ylim([1.5 2])
xlim([1 100])
xticks_init = xticks; 
xticks([1, xticks_init])
exportgraphics(ax, 'Variable_stepsizes_k=-1e-03.pdf')
%% Compare denomniators (semilogy)
figure(2); clf; leg_names = {};
semilogy(Nvec, 1./P_den_hs, 'm', 'LineWidth', 4); hold on; leg_names{end+1} = 'Asymptotically $(\gamma L)_*(\kappa)$';
semilogy(Nvec, 1./P_den_hs_N,'-b', 'LineWidth',4); leg_names{end+1} = "Best constant $(\gamma L)_*(N,\kappa)$";
semilogy(Nvec, 1./P_den_seq,'k', 'LineWidth',4); leg_names{end+1} = "$s_{N}(\kappa)$";
semilogy(Nvec, 1./P_den_h_dyn, 'r:', 'LineWidth', 4); hold on; leg_names{end+1} = 'Scheduled';
legend(leg_names, 'Interpreter','latex','FontSize',17, 'Location','northeast')
xlabel('Number of iterations $N$')
ylabel('Worst-case guarantee', 'Interpreter','latex')
title_str = sprintf("Denominators of rates ($\\kappa=%.e$)", kappa);
title(title_str, 'Interpreter','latex')
xlim([0 300])
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 17;
% ax.GridLineWidth = 1;
grid on
exportgraphics(ax, 'Dens_Variable_stepsizes_k=-1e-03.pdf')
% return
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