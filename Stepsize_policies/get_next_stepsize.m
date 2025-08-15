function hp = get_next_stepsize (h, k)
    if nargin == 0
        h = 0; k = .1;
    end
    hp_candidates = roots([-k*(2-h*(1+k)), 2*(1+k)*(2-h*(1+k)), 2*(-3+2*h*(1+k)), -2*h]);
    hp = hp_candidates ( (hp_candidates <= 2/(1+k) + 10*eps ) & hp_candidates > 0 );    
    if isempty(hp)
        disp("Problem")
    end
end