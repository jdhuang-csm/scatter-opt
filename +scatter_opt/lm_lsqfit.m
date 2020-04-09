function [fit,internal] = lm_lsqfit(data,L_air,L,mu_x0,eps_x0,wg_x0,...
    muNPoles1,epsNPoles1,NPoles2,varargin)
% Parameters
% ----------
% data: table with frequency and S_ij
% L_tot : total sample holder thickness
% mu_x0 : initial guess parameter vector for Laurent model for mu
% eps_x0 : initial guess parameter vector for Laurent model for epsilon
% wg_x0 : initial guess parameters for waveguide. Entries are L1, L2, and
%   lambda_c
% Named (optional) parameters
% ----------------
% MaxFEval : maximum number of function evaluations to allow the optimizer
%   to perform
% MaxIter : maximum number of iterations to allow the optimizer to perform
% MultiStart: if true, start the optimization from multiple initial
%   points
% MultiStartPoints : number of points from which to initialize the
%   optimization
% ShiftFOPoles : if true, shift the parameters for first-order poles when 
%   generating MultiStart points. If false, leave the first-order pole
%   parameters the same for all MultiStart points (only shift 2nd-order
%   poles)
% NCCP2 : number of complex conjugate pairs in second-order poles. 0 (real 
%   2nd-order poles only) tends to work well, but in some cases complex
%   2nd-order poles may be beneficial.
    import scatter_opt.*
	parser = inputParser;
    addParameter(parser,'MaxFEval',-1)
    addParameter(parser,'MaxIter',-1)
    addParameter(parser,'MultiStart',true)
    addParameter(parser,'MultiStartPoints',10)
    addParameter(parser,'ShiftFOPoles',false)
    addParameter(parser,'NCCP2',0)
    addParameter(parser,'UseSP',{'s11' 's21' 's12' 's22'})
    parse(parser,varargin{:})
    
    % lsqcurvefit cannot handle complex functions with bounds. Thus, to 
    % simultaneously optimize s11, s21, s12, and s22, need to split each
    % s_ij into its real and imaginary parts and concatenate all
    % components into a vector of length 8*N. Similarly, the input  
    % frequencies must be duplicated 8 times to match the s_ij vector.
    use_sp = parser.Results.UseSP;
    xdata = repmat(data.freq,length(use_sp)*2,1);
    ydata = zeros(length(xdata),1);
    for j=1:length(use_sp)
        sp = use_sp{j};
        nfreq = length(data.freq);
        ydata((j-1)*2*nfreq+1:j*2*nfreq) = [real(data.(sp)); imag(data.(sp))];
    end
%     ydata = [real(data.s11); imag(data.s11);...
%         real(data.s21); imag(data.s21);...
%         real(data.s12); imag(data.s12);...
%         real(data.s22); imag(data.s22)];
    
    % Prepare mu_x0 and eps_x0 for lsq fit
    % Originally, was going to test different numbers of complex conjugate 
    % pairs in 2nd-order poles, but it seems that real 2nd-order
    % poles give similar or better solutions and converge faster (not
    % always though...)
    num_ccp = parser.Results.NCCP2;
    if num_ccp > NPoles2/2
        error('Number of 2nd-order complex conjugate pairs cannot exceed NPoles/2')
    end
    
    if isempty(mu_x0)
        mu_x0_prep = [];
        mu_np1 = 0;
        mu_np2 = 0;
    else
        [mu_x0_prep,mu_np1,mu_np2] = lm_prep_x0(mu_x0,num_ccp,muNPoles1,NPoles2);
    end
    [eps_x0_prep,eps_np1,eps_np2] = lm_prep_x0(eps_x0,num_ccp,epsNPoles1,NPoles2);
    
    % Need to simultaneously optimize parameter vectors for mu and eps.
    % Create a single parameter vector with mu, eps, and waveguide 
    % parameters
    lenx_mu = length(mu_x0_prep);
    lenx_eps = length(eps_x0_prep);
    x0 = [mu_x0_prep eps_x0_prep wg_x0];
    
    % get bounds for mu_x and eps_x
    if mu_np1==0
        lb_xmu = []; ub_xmu = [];
    else
        [lb_xmu,ub_xmu] = lm_xbounds(mu_x0_prep,mu_np1,mu_np2);
    end
    [lb_xeps,ub_xeps] = lm_xbounds(eps_x0_prep,eps_np1,eps_np2);
    % set bounds for waveguide parameters
    % L should be within 15% of estimate
    lb_L = L*0.85; ub_L = L*1.15;
    % L1 should be <= 15% of L
    lb_L1 = -L*0.15; ub_L1 = L*0.15;
    % lambda_c should be +/- 10% of estimate
    lambda_c = wg_x0(3);
    lb_lamc = lambda_c*0.9; ub_lamc = lambda_c*1.1;
    % concatenate waveguide bounds
    lb_wg = [lb_L lb_L1 lb_lamc];
    ub_wg = [ub_L ub_L1 ub_lamc];
    % concatenate all bounds
    x_lb = [lb_xmu lb_xeps lb_wg];
    x_ub = [ub_xmu ub_xeps ub_wg];

    % prepare for optimization
    % create anonymous function to pass args to lm_lsqfun
    fun = @(x,xdata)lm_lsqfun(x,xdata,L_air,lenx_mu,lenx_eps,...
        mu_np1,mu_np2,eps_np1,eps_np2,parser.Results.UseSP);
    % set optimization options
    typical_vals = x0;
    % set L and L1 typical values
    typical_vals(end-2) = L;
    typical_vals(end-1) = L_air*0.05;
    % a0 and imag components for real-valued a2 and b2 will still be zero. 
    % Set typical vals for these to 1
    typical_vals(typical_vals==0)=1;
    options = optimoptions('lsqcurvefit','TypicalX',typical_vals);
    if parser.Results.MaxFEval~=-1
        options.MaxFunctionEvaluations = parser.Results.MaxFEval;
    end
    if parser.Results.MaxIter~=-1
        options.MaxIterations = parser.Results.MaxIter;
    end
    
    % optimize x
    if parser.Results.MultiStart
        rng default
        problem = createOptimProblem('lsqcurvefit','objective',fun,...
            'x0',x0,'xdata',xdata,'ydata',ydata,'lb',x_lb,'ub',x_ub,...
            'options',options);
        ms = MultiStart('StartPointsToRun','bounds',...
            'XTolerance',1e-3,'Display','iter');
        
        % generate start points
        startpts = repmat(x0,parser.Results.MultiStartPoints,1);
        
        % randomly shift start values of 2nd-order poles by up to a factor
        % of 10
        mu2_index = (3+mu_np1*4:2+(mu_np1+mu_np2)*4);
        eps2_index = (3+lenx_mu+eps_np1*4:2+lenx_mu+(eps_np1+eps_np2)*4);
        mu2_points = lm_mspoints(x0(mu2_index),...
            parser.Results.MultiStartPoints-1,10);
        eps2_points = lm_mspoints(x0(eps2_index),...
            parser.Results.MultiStartPoints-1,10);
        startpts(2:end,mu2_index) = mu2_points;
        startpts(2:end,eps2_index) = eps2_points;
        
        % Careful with shifting 1st-order poles - even tiny shifts (<1%)
        % can lead to much higher residuals in the solution. Shifting
        % 1st-order poles tends to be helpful when optimizing a larger
        % number of poles (i.e. 3+ first-order poles)
        if parser.Results.ShiftFOPoles
            % apply smaller shifts to 1st-order poles (up to factor of 1.1)
            mu1_index = (1:2+mu_np1*4);
            eps1_index = (lenx_mu+1:lenx_mu+2+eps_np1*4);
            mu1_points = lm_mspoints(x0(mu1_index),...
                parser.Results.MultiStartPoints-1,1.2);
            eps1_points = lm_mspoints(x0(eps1_index),...
                parser.Results.MultiStartPoints-1,1.2);
            startpts(2:end,mu1_index) = mu1_points;
            startpts(2:end,eps1_index) = eps1_points;
        end
        
        startpts = CustomStartPointSet(startpts);
        [x_opt,resnorm,exitflag,output,solutions] = run(ms,problem,...
            startpts);
    else
        [x_opt,resnorm,exitflag,output] = lsqcurvefit(fun,x0,xdata,...
            ydata,x_lb,x_ub,options);
    end
    
    % separate x for mu, eps, and wg
    x_mu = x_opt(1:lenx_mu);
    x_eps = x_opt(lenx_mu+1:lenx_mu+lenx_eps);
    x_wg = x_opt(lenx_mu+lenx_eps+1:end);
    % reconstruct x_mu and x_eps with all complex conjugate pairs
    if mu_np1==0
        x_mu = [];
        mnp1 = 0; mnp2 = 0;
    else
        [x_mu,mnp1,mnp2] = lm_repair_x(x_mu,mu_np1,mu_np2);
    end
    [x_eps,enp1,enp2] = lm_repair_x(x_eps,eps_np1,eps_np2);
    
    % put results in fit struct
    fit = struct;
    fit.x_mu = x_mu;
    fit.mu_np1 = mnp1;
    fit.mu_np2 = mnp2;
    fit.x_eps = x_eps;
    fit.eps_np1 = enp1;
    fit.eps_np2 = enp2;
    fit.L = x_wg(1);
    fit.L1 = x_wg(2);
    fit.L2 = L_air - (x_wg(1) + x_wg(2));
    fit.lambda_c = x_wg(3);
    fit.resnorm = resnorm;
    
    % store internal data for troubleshooting
    internal = struct;
    internal.x_opt = x_opt;
    internal.xdata = xdata;
    internal.ydata = ydata;
    internal.lenx_mu = lenx_mu;
    internal.lenx_eps = lenx_eps;
    internal.x_lb = x_lb;
    internal.x_ub = x_ub;
    internal.mu_np1 = mu_np1;
    internal.mu_np2 = mu_np2;
    internal.eps_np1 = eps_np1;
    internal.eps_np2 = eps_np2;
    internal.lsq_exitflag = exitflag;
    internal.lsq_output = output;
    if parser.Results.MultiStart
        internal.lsq_solutions = solutions;
    end
end
    
    
        
    
    
    
    
    