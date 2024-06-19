clear; clc; close all;
% Script to plot the denominator in the rates with constant stepsizes for
% full range (0,2), with fixed:
% - number of iterations: N
% - initial function gap: Delta
% - upper curvature: L=1
% - curvautre ratios: kappa_vec
%%
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
L = 1; 
N = 10; Delta = 1;
kappa_vec = [-1, -0.5, -0.1, -0.01, -0.005, -0.001, 0];
figure(1); clf; leg_names = {};
for kk = 1 : length(kappa_vec)
    kappa = kappa_vec(kk);
    mu = kappa*L;
    % compute first N thresholds
    gL_thr_vec = zeros(1,N);
    Nvec = 1 : N;
    for gg = 1 : length(gL_thr_vec)
        gL_thr_vec(gg) = get_gammaL_bar(kappa,Nvec(gg));
    end    
    % define the stepsizes vector 
    gL_vec = sort(unique([gL_thr_vec, linspace(0,gL_thr_vec(1),100), linspace(gL_thr_vec(1),gL_thr_vec(2),1000), linspace(gL_thr_vec(2),2,5000)]));
    wc_rate = 0 * gL_vec;
    DS = 0 * gL_vec;
    for gg = 1 : length(gL_vec)
        [wc_rate(gg), ~, DS(gg), ~] = get_wc_rate(L,mu,gL_vec(gg),N,Delta);   
    end
    plot(gL_vec, DS, 'LineWidth', 2); hold on; grid on
    leg_names{end+1} = sprintf('$\\kappa = %.3f$', kappa);
end
legend(leg_names, 'Interpreter','latex', 'Location','northwest');
ylabel("Denominator $\frac{2L[f(x_0)-f(x_N)]}{\min \{\nabla f(x_i)\|^2\}}$");
xlabel("Constant stepsizes $\gamma L$");
ax = gca; ax.FontWeight = 'bold'; ax.FontSize = 16; ax.GridLineWidth = 2;