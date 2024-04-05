function hi_vec = get_sequence_stepsizes (k, N)
    if nargin == 0
        k = .1; N = 5;
    end
    hi_vec = zeros(1, N+1);
    for i = 1 : N
        hi_vec(i+1) = get_next_stepsize (hi_vec(i), k);
    end
end