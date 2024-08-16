# Numerics accompanying tight rates for gradient descent
 ## Teodor Rotaru, François Glineur, and Panagiotis Patrinos. Exact worst-case convergence rates of gradient descent: a complete analysis for all constant stepsizes over nonconvex and convex functions, 2024; URL https://arxiv.org/abs/2406.17506.
 
Folders and files:

- [main folder]
    - gradient_descent_PEP.m: Solves the PEP on Gradient Descent for Smooth functions: 
	    - To run it, PESTO toolbox must be installed and updated with the files from folder PESTO_files (can be copy-pasted directly over folder PESTO_files from PESTO toolbox):
		    - pep.m -- adapted to include hypoconvex functions;
		    - Functions_classes -> Hypoconvex.m -- implements the class of hypoconvex functions.
        - It confirms numerically the tightness for all the rates (Theorems 2.1-2.4).
    - GD_rates_live_script.mlx: Matlab live script (~notebook) to check the algebraic manipulations for the upper bounds proofs

- Tightness_check_triplets:
        - build_2D_interpolating_triplets_variable.m: Tightness of Conjecture 1 (variable step-sizes in [1, \bar{\gamma}^1];
        - build_3D_interpolating_triplets.m: Tightness of Conjecture 2 (fixed step-sizes in [\bar{\gamma}^1, \bar{\gamma}^{N-1}];

- Compute_rates: it includes functions to compute the theoretical worst-case rates from Theorems 2.1-2.4
    - get_wc_rate.m: computes the rates with:
        - constant stepsizes (Theorems 2.1-2.3);
        - variable stepsizes up to threshold (\bar{\gamma}L)^1 (Theorem 2.4).
    - get_gammaL_bar.m: compute the stepsize thresholds \bar{\gamma}^i; (see Section 2.5)
    - get_Nbar.m: for fixed stepsize gamma, it computes the maximum number of iterations for which the worst-case rate is (1-\gamma L)^{2N}, corresponding to worst-case function Lx^2/2;
    - get_PN.m: it computes the (complicated) denominator for hypoconvex functions.
	- plot_denominators_for_constant_stepsizes.m: it shows the complicated denominator in the rate for hypoconvex functions and constant stepsize belonging to (0, 2).
		
- Stepsize_policies (Section 2.6)
    - it includes files to compute the optimal constant stepsizes (Section 2.6.1) and dynamic stepsize sequences (Section 2.6.2)
      - compute_optimal_stepsize_hypoconvex.m: (Section 2.6.1, Proposition 2.15 and used to generate Figure 2.5)
            - computes the optimal constant stepsize (possibly dependent on the number of iterations) for hypoconvex (weakly convex) functions (used for comparisons);
      - get_sequence_stepsizes.m: (Definition 2.16)
            - generates the dynamic sequence (independent on the number of iterations) converging to 2/(1+max{0,kappa}) for nonconvex and strongly convex functions;
      - variable_stepsizes_comparison.m: (used to generate Tables 2.1 and 2.2)
            - compares, for convex / strongly convex functions, the rates with different stepsize schedules:
      - variable_stepsizes_comparison_hypoconvex.m: compares, for hypoconvex (weakly convex) functions, the rates with different stepsize schedules;
      - gradient_descent_dynamic_stepsizes_PEP.m: (confirms numerically Theorem 2.18, Corollary 2.19, Theorem 2.20, Proposition 2.21)
            - solves the PEP for Gradient Descent using dynamic stepsizes; 
            - for nonconvex functions, the sequence is trimmed by the (asymptotic) optimal constant stepsize;
            - confirms the tightness of rates obtained using dynamic stepsizes.
