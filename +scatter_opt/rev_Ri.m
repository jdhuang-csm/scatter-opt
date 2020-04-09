function R = rev_Ri(omega,lambda_c,Li)
% rotation term
    import scatter_opt.*
	R = exp(-rev_gamma_0(omega,lambda_c).*Li);
end