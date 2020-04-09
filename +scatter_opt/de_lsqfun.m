function F = de_lsqfun(x,xdata,L_air)
    import scatter_opt.*
	x_eps = x(1:end-3);
%     [p_vec,z_vec,C] = de_expand_x(x_eps);
    x_wg = x(end-2:end);
    L = x_wg(1);
    L1 = x_wg(2);
    lambda_c = x_wg(3);
    L2 = L_air - (L + L1);
%     disp(L);disp(L1);disp(L2);disp(lambda_c)

    % get predicted epsilon
    freq = xdata(1:round(length(xdata)/8));
    eps = de_eval(x_eps,freq);
    % mu set to unity
    mu = ones(length(freq),1);
    
    % perform the reverse transformation to get s_ij
    rev_in = array2table([freq mu eps],'VariableNames',{'freq' 'mu' 'eps'});
    rev_out = rev_transform(rev_in,lambda_c,L,L1,L2);

    % concatenate the s_ij vectors
    F = [real(rev_out.s11); imag(rev_out.s11);...
        real(rev_out.s21); imag(rev_out.s21);...
        real(rev_out.s12); imag(rev_out.s12);...
        real(rev_out.s22); imag(rev_out.s22)];
end