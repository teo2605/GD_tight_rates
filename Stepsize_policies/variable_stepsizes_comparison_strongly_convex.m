clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% setup
kappa = 1e-4; L = 1; Delta = 1; % Delta = f0 - f*
Nvec = [1,2,5,10,20,30,40,50,70,100];
% Nvec = 1 : 100;
Nvec_sel = Nvec + 1;
Nmax = max(Nvec);
%% get stepsizes
%% variable stepsizes
h_var_vec = get_sequence_stepsizes (kappa, Nmax);
%% optimal constant stepsize = the N-th stepsize threshold
h_ct_opt_vec = zeros(1, Nmax+1);
for N = 1 : Nmax
    h_ct_opt_vec(N+1) = get_gammaL_bar(kappa,N);
end
%% rates
wc_var_stepsizes = @(h) ( (2-h*(1+kappa)) ./ (2-h*kappa) );
wc_var_stepsizes_N_its = wc_var_stepsizes(h_var_vec);
wc_var_stepsizes_N_its_inv = 1./wc_var_stepsizes(h_var_vec);
wc_ct_stepsizes = (1-h_ct_opt_vec).^(2*(0:Nmax));
h_celebrated = 2/(1+kappa);
wc_celebrated = (1-h_celebrated).^(2*(0:Nmax));
h_celebrated_1oL = 1;
wc_celebrated_1oL = 1./(1 + h_celebrated_1oL * ( (-1+(1-kappa*h_celebrated_1oL).^(-2*(0:Nmax))) / (h_celebrated_1oL*kappa) ) );
%% plot variable step-size
figure(1); clf; leg_names={};
plot(0:Nmax, 2/(1+kappa)*ones(1,Nmax+1), '-k', 'LineWidth', 4); leg_names{end+1} = 'Celebrated $\frac{2}{1+\kappa}$';
hold on;  grid on;
plot(0:Nmax, h_ct_opt_vec, '-b', 'LineWidth', 4); leg_names{end+1} = 'Best constant $\overline{\gamma L}_N(\kappa)$';
plot(0:Nmax, h_var_vec, ':r', 'LineWidth', 4); leg_names{end+1} = 'Scheduled $s_N(\kappa)$';

title_str = sprintf("Sequence of stepsizes for $\\kappa=%.e$", kappa);
% title(title_str, 'Interpreter','latex')
legend(leg_names, 'Interpreter','latex', 'Location', 'southeast', 'fontsize', 17)
xlabel('Number of iterations $N$', 'Interpreter','latex')
ylabel('Stepsize $\gamma L$', 'Interpreter','latex')
xlim([1 100]); %ylim([1.5 2])
ylim([1.48, 2])
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 17;
% ax.GridLineWidth = 2;

xticks_init = xticks; 
xticks([1, xticks_init])
exportgraphics(ax, 'Variable_stepsizes_k=1e-04.pdf')
%% plot rates
figure(2); clf; leg_names={};
  
semilogy(0:Nmax, wc_celebrated, '-k', 'LineWidth', 4); leg_names{end+1} = 'Celebrated $\frac{2}{1+\kappa}$'; grid on; hold on;
semilogy(0:Nmax, wc_celebrated_1oL, '-m', 'LineWidth', 4); leg_names{end+1} = 'Celebrated $1$';
semilogy(0:Nmax, wc_ct_stepsizes, '-b', 'LineWidth', 4); leg_names{end+1} = 'Best constant $\overline{\gamma L}_N(\kappa)$';
semilogy(0:Nmax, wc_var_stepsizes_N_its, ':r', 'LineWidth', 4); leg_names{end+1} = 'Scheduled $s_N(\kappa)$';

title_str = sprintf("Convergence rates ($\\kappa=%.e$)", kappa);
% title(title_str, 'Interpreter','latex')
legend(leg_names, 'Interpreter','latex', 'Location', 'northeast', 'fontsize', 17)
xlabel('Number of iterations $N$', 'Interpreter','latex')
ylabel('Worst-case guarantee', 'Interpreter','latex')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 17;
% ax.GridLineWidth = 1;
exportgraphics(ax, 'Dens_Variable_stepsizes_k=1e-04.pdf')
% return
%% plot rates' denominators
figure(3); leg_names={};
semilogy(0:Nmax, 1./wc_var_stepsizes_N_its-1, '--ro', 'LineWidth', 2); grid on; hold on; leg_names{end+1} = 'Scheduled $s_N(\kappa)$';
semilogy(0:Nmax, 1./wc_celebrated-1, '-k*', 'LineWidth', 2); leg_names{end+1} = 'Celebrated $\frac{2}{1+\kappa}$';
semilogy(0:Nmax, 1./wc_ct_stepsizes-1, '-b*', 'LineWidth', 2); leg_names{end+1} = 'Best constant $\overline{\gamma L}_N(\kappa)$';
title_str = sprintf("Denominators of rates ($\\kappa=%.e$)", kappa);
title(title_str, 'Interpreter','latex')
legend(leg_names, 'Interpreter','latex', 'Location', 'southeast')
xlabel('Number of iterations $N$', 'Interpreter','latex')
ylabel('Denominators', 'Interpreter','latex')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
% ax.GridLineWidth = 2;
%% Table: Comparison of denominators
Iterations = num2str(Nvec');
var_stepsizes = h_var_vec(Nvec_sel)';
Den_var_steps = (1./wc_var_stepsizes_N_its(Nvec_sel))';
ct_stepsizes = h_ct_opt_vec(Nvec_sel)';
Den_ct_steps = (1./wc_ct_stepsizes(Nvec_sel))';
Ratio = Den_var_steps ./ Den_ct_steps * 100;
wc_rate_var_steps = 2*L*Delta*wc_var_stepsizes_N_its(Nvec_sel)';
rate_comparison = table( Iterations, ct_stepsizes, Den_ct_steps, var_stepsizes, Den_var_steps, Ratio, wc_rate_var_steps);
n_decimal = 5;
table_compare_rates = varfun(@(x) num2str(x, ['%' sprintf('.%df', n_decimal)]), rate_comparison);
table_compare_rates.Properties.VariableNames = rate_comparison.Properties.VariableNames;
table_compare_rates.Properties.RowNames = rate_comparison.Properties.RowNames;
disp(table_compare_rates)