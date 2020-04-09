function out = de_rev_eval(defit,freq)
	import scatter_opt.*
    mu = ones(length(freq),1);  
    eps = de_eval(defit.x_eps,freq);
    
    % perform the reverse transformation to get s_ij
    rev_in = array2table([freq mu eps],'VariableNames',...
        {'freq' 'mu' 'eps'});
    out = rev_transform(rev_in,defit.lambda_c,defit.L,defit.L1,defit.L2);
end