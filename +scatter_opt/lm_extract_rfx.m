function x = lm_extract_rfx(fit)
% Extract parameter vector from rationalfit object
    import scatter_opt.*
	a0 = fit.D;
    a1 = fit.C;
    b1 = fit.A;
    num_poles = length(a1);
    a2 = zeros(1,num_poles);
    b2 = zeros(1,num_poles);
    x = lm_flatten_terms(a0,a1,b1,a2,b2);
end