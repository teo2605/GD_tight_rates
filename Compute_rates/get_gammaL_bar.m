function [gL_thr, it, tol] = get_gammaL_bar(k, i, gL0, tol)
%% 
% input: k = mu/L; i - integer threshold
% output: gL_thr - limit step-size between transition to next i
% Observation: Implemenation via standard Newton method
%%
if nargin == 0
    % k = 0.743774; 	 i = 6 ;
    k = -1; i = 10; 
end
if nargin < 3        
    gL0 =  3 / ( 1 + k + sqrt(1-k+k^2) ); % initial condition corresponding
end
if nargin < 4
    tol = 1e-9; 
    tol = eps;
end
if i == 0;  gL_thr = 1; it = 0;  return;  end
if i == 1;  gL_thr = 3 / (1+k+sqrt(1-k+k^2)); it = 0; return; end
%% implement Newton method to find the root of Ti
itmax = 1000; 
f  = @(gL) get_Ti  (gL,k*gL,i);
df = @(gL) get_d_Ti(gL,k*gL,i);

x_max = 2/(1+max(0,k));
x_hi = x_max;
x_lo = gL0;
x = ( x_lo + x_hi )/2; it = 0;
while ~( (abs(f(x)) < tol) ) && it < itmax && ~( (df(x))^2 < tol) && (x_hi - x_lo > eps)
    if f(x) > tol % f(x) > 0
        x_hi = x;
    else % f(x) < 0
        x_lo = x;
    end    
    x = ( x_lo + x_hi )/2;
    it = it+1;
    % [it, x, f(x), x_lo, f(x_lo), x_hi, f(x_hi)]
end
if it == itmax 
    % fprintf("x0 = %.4f \t x = %.4f \t f(x) = %.4e \t |df(x)| = %.4e \n", gL0, x, f(x), df(x));
    fprintf("x0 = %.6f \t x = %.6f \t k = %.6f \t Nbar = %d \t f(x) = %.6e  \t |df(x)| = %.4e \n", gL0, x, k, i, f(x), df(x));
    error("Exceeded number of iterations!");
end
gL_thr = x;
end
%% auxiliary functions
function dTi_dgL = get_d_Ti(gL,gmu,i)   
    dTi_dgL = -get_Ti(gL,gmu,i)/gL + 2*i/gL * ( (1-gmu)^(-2*i-1) - (1-gL)^(-2*i-1)  );
end