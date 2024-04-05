function [Nbar, it] = get_Nbar(gL,gmu,tol)    
    if nargin == 0
        gL = 1.75;
        k = .1;
        gmu = gL * k;
    end
    if nargin < 3; tol = 1e-12; end
    k = gmu / gL;
    if  gL < 3 / ( 1 + k + sqrt(1-k+k^2) )
        Nbar = 0; return;
    end
    n = 1-gL;    m = 1-gmu;
    if m+n <= 0 
        Nbar = Inf; 
        return
    end
    %% if gmu < 0 (hypoconvex)
    Nbar_lb_n = log( 1 + (1-n)/(m+n) ) / ( log (1/n^2) );  % lower bound wrt. gL
    Nbar_lb_m = log( 1 + (1-m)/(m+n) ) / ( log (1/m^2) );  % lower bound wrt. gmu
    if ~isreal ( [Nbar_lb_n, Nbar_lb_m] )
        Nbar_lb = 0; 
    else
        Nbar_lb = ceil( -1 + max ( Nbar_lb_m, Nbar_lb_n ));
    end
    check_Nbar        = @(Nbar) (get_Ti(gL,gmu,Nbar  ) >= -tol); % check T(gL,gmu,Nbar  ) >=0
    check_Nbar_p      = @(Nbar) (get_Ti(gL,gmu,Nbar+1) < tol); % check T(gL,gmu,Nbar+1) < 0
    check_Nbar_Nbar_p = @(Nbar) (check_Nbar(Nbar)  && check_Nbar_p (Nbar)); % check that both conditions are satisfied
    %% Apply binary search
    itmax = 1e6; it = 1;
    Nbar = Nbar_lb; % initial guess
    if ~check_Nbar_Nbar_p(Nbar_lb)        
        Nbar_ub = Nbar_lb + 1;
        it_bnd = 0;
        %% find feasible upper bound 
        while ~check_Nbar_p(Nbar_ub)
            Nbar_lb = Nbar_ub;
            Nbar_ub = 2*Nbar_ub;
            it_bnd = it_bnd + 1;
            if it_bnd > itmax
                fprintf("%.4f, %.4f, %.4f\n", gL, gmu, Nbar)
                Nbar = inf;
                return;
                % error("Exceeded number of iteration to compute initial bounds of Nbar!");
            end
        end
        %% find Nbar by binary search
        while ~check_Nbar_Nbar_p(Nbar)
            if check_Nbar_Nbar_p(Nbar_ub)
                Nbar = Nbar_ub; return;
            end
            if check_Nbar_Nbar_p(Nbar_lb)
                Nbar = Nbar_lb; return;
            end
            Nbar = floor( (Nbar_lb + Nbar_ub)/2 );
            if ~check_Nbar_p(Nbar) % Nbar too low
               Nbar_lb = Nbar;               
            else
                if ~check_Nbar (Nbar) % Nbar too high
                    Nbar_ub = Nbar;
                else
                    break;
                end
            end
            it = it + 1;
            if it > itmax
                fprintf("%.4f, %.4f, %.4f\n", gL, gmu, Nbar)
                error("Exceeded number of iteration to compute Nbar!");
            end
        end
    end
end