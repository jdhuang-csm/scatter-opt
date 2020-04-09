function out = lm_rev_eval(lmfit,freq)
% Evaluate s-parameters based on lmfit
% Parameters
% ----------
% lmfit: lmfit object from lm_lsqfit
% freq: frequencies to evaluate
    import scatter_opt.*
	mu = lm_eval(lmfit.x_mu,freq,lmfit.mu_np1,lmfit.mu_np2);
    eps = lm_eval(lmfit.x_eps,freq,lmfit.eps_np1,lmfit.eps_np2);
%     disp(lmfit.mu_np1);disp(lmfit.mu_np2)
%     disp(lmfit.eps_np1);disp(lmfit.eps_np2)
%     disp(freq(1:10))
%     disp(mu(1:10))
%     disp(eps(1:10))
    rev_in = array2table([freq mu eps],'VariableNames',{'freq' 'mu' 'eps'});

    % perform the reverse transformation to get s_ij
    out = rev_transform(rev_in,lmfit.lambda_c,lmfit.L,lmfit.L1,lmfit.L2);
end
    