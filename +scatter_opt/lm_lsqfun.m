function F = lm_lsqfun(x,xdata,L_air,lenx_mu,lenx_eps,mu_np1,mu_np2,...
    eps_np1,eps_np2,use_sp)
% Function for use in lsqcurvefit. Evaluate s11, s21, s12, and s22 based on
% parameter vector x and input frequencies xdata.
    import scatter_opt.*
	% extract mu and eps parameter vectors
    if lenx_mu > 0
        x_mu = x(1:lenx_mu);
        [x_mu,mnp1,mnp2] = lm_repair_x(x_mu,mu_np1,mu_np2);
    end
    x_eps = x(lenx_mu+1:lenx_mu+lenx_eps);
    % re-pair unpaired complex conjugates
%     [a0,a1,b1,a2,b2] = lm_expand_x(x_mu,mu_np1,mu_np2)
    
    [x_eps,enp1,enp2] = lm_repair_x(x_eps,eps_np1,eps_np2);
%     disp(x_mu)
%     disp(x_eps)
%     disp('x_mu')
%     fprintf([repmat('%7.4e    ',1,size(x_mu,2)),'\n'],x_mu')
%     disp('x_eps')
%     fprintf([repmat('%7.4e    ',1,size(x_eps,2)),'\n'],x_eps')
    % extract waveguide parameters
    x_wg = x(lenx_mu+lenx_eps+1:end);
    L = x_wg(1);
    L1 = x_wg(2);
    lambda_c = x_wg(3);
    L2 = L_air - (L + L1);
%     disp(L);disp(L1);disp(L2);disp(lambda_c)

    % get predicted mu and epsilon
    freq = xdata(1:round(length(xdata)/8));
%     disp(mu_np1);disp(mu_np2)
%     disp(eps_np1);disp(eps_np2)
    if lenx_mu==0
        mu = ones(length(freq),1);
    else
        mu = lm_eval(x_mu,freq,mnp1,mnp2);
    end
    eps = lm_eval(x_eps,freq,enp1,enp2);
%     disp(freq(1:10))
%     disp(mu(1:10))
%     disp(eps(1:10))
    rev_in = array2table([freq mu eps],'VariableNames',{'freq' 'mu' 'eps'});

    % perform the reverse transformation to get s_ij
    rev_out = rev_transform(rev_in,lambda_c,L,L1,L2);

    % concatenate the s_ij vectors
    F = zeros(length(xdata),1);
    for j=1:length(use_sp)
        sp = use_sp{j};
        nfreq = length(freq);
        F((j-1)*2*nfreq+1:j*2*nfreq) = [real(rev_out.(sp)); imag(rev_out.(sp))];
    end
%     F = [real(rev_out.s11); imag(rev_out.s11);...
%         real(rev_out.s21); imag(rev_out.s21);...
%         real(rev_out.s12); imag(rev_out.s12);...
%         real(rev_out.s22); imag(rev_out.s22)];
%     disp(length(F))
end