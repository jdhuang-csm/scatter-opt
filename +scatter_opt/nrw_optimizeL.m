function [L1,L2] = nrw_optimizeL(Sdata,L,L_air,lambda_c)
% Determine values of L1 and L2 that most closely match s11 and s22 across
% the measurement frequency range
    import scatter_opt.*
    gam0 = rev_gamma_0(Sdata.freq*2*pi,lambda_c);
    function C = objfun(x,Sdata,gam0,L,L_air,reg_coef)
        L1_ = x(1); L2_ = x(2); 
        L_ = L_air - L1_ - L2_;
        R1 = exp(-gam0*L1_);
        R2 = exp(-gam0*L2_);
        s11_s = Sdata.s11./(R1.^2);
        s22_s = Sdata.s22./(R2.^2);
        C = sum(abs(s11_s-s22_s).^2);
        % penalty for departing from measured value of L
        penalty = reg_coef*(L_-L).^2;
        C = C + penalty;
    end
    
    x0 = [(L_air-L)/2 (L_air-L)/2];
    fun = @(x)objfun(x,Sdata,gam0,L,L_air,1e3);
    [x_opt,rss] = fminsearch(fun,x0);
    
    L1 = x_opt(1); L2 = x_opt(2);
end
    