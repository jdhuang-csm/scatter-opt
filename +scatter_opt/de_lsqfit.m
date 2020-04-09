function [fit,internal] = de_lsqfit(data,L_air,L,eps_x0,wg_x0,varargin)
% Parameters
% ----------
% data: table with frequency and S_ij
% L_tot : total sample holder thickness
% eps_x0 : initial guess parameter vector for epsilon model
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
    import scatter_opt.*
	parser = inputParser;
    addParameter(parser,'MaxFEval',-1)
    addParameter(parser,'MaxIter',-1)
    addParameter(parser,'MultiStart',true)
    addParameter(parser,'MultiStartPoints',10)
    parse(parser,varargin{:})
    
    % lsqcurvefit cannot handle complex functions with bounds. Thus, to 
    % simultaneously optimize s11, s21, s12, and s22, need to split each
    % s_ij into its real and imaginary parts and concatenate all
    % components into a vector of length 8*N. Similarly, the input  
    % frequencies must be duplicated 8 times to match the s_ij vector.
    xdata = repmat(data.freq,8,1);
    ydata = [real(data.s11); imag(data.s11);...
        real(data.s21); imag(data.s21);...
        real(data.s12); imag(data.s12);...
        real(data.s22); imag(data.s22)];
    
    % Create a single parameter vector with eps and waveguide parameters
    x0 = [eps_x0 wg_x0];
    
    % get bounds for mu_x and eps_x
    [lb_xeps,ub_xeps] = de_xbounds(eps_x0);
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
    x_lb = [lb_xeps lb_wg];
    x_ub = [ub_xeps ub_wg];
    
    % prepare for optimization
    % create anonymous function to pass args to lm_lsqfun
    fun = @(x,xdata)de_lsqfun(x,xdata,L_air);
    % set optimization options
    typical_vals = x0;
    % set L and L1 typical values
    typical_vals(end-2) = L;
    typical_vals(end-1) = L_air*0.05;
    % shouldn't be any zeros remaining in x0 vector
%     typical_vals(typical_vals==0)=1;
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
        
        % randomly shift start values of zeros, poles, and C by up to a 
        % factor of 10
%         eps_index = (1:length(eps_x0));
        n_poles = (length(eps_x0)-2)/3;
        p_index = (1:n_poles);
        zc_index = (1+n_poles:length(eps_x0));
        % allow poles to shift by factor of 2
        p_points = lm_mspoints(eps_x0(p_index),...
            parser.Results.MultiStartPoints-1,4);
        % allow other params to shift by factor of 10
        zc_points = lm_mspoints(eps_x0(zc_index),...
            parser.Results.MultiStartPoints-1,10);
        startpts(2:end,p_index) = p_points;
        startpts(2:end,zc_index) = zc_points;
        
        startpts = CustomStartPointSet(startpts);
        [x_opt,resnorm,exitflag,output,solutions] = run(ms,problem,...
            startpts);
    else
        [x_opt,resnorm,exitflag,output] = lsqcurvefit(fun,x0,xdata,...
            ydata,x_lb,x_ub,options);
    end
    
    % separate x_opt for eps and wg
    x_eps = x_opt(1:length(eps_x0));
    x_wg = x_opt(1+length(eps_x0):end);
    
    % put results in fit struct
    fit = struct;
    fit.x_eps = x_eps;
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
    internal.x_lb = x_lb;
    internal.x_ub = x_ub;
    internal.lsq_exitflag = exitflag;
    internal.lsq_output = output;
    if parser.Results.MultiStart
        internal.lsq_solutions = solutions;
    end
end