function [mu_a2,mu_b2,eps_a2,eps_b2] = lm_initialfit2(data,mu_x01,eps_x01,varargin)
% Estimate terms for 2nd-order poles    
	import scatter_opt.*
	parser = inputParser;
    addParameter(parser,'NPoles',1)
    addParameter(parser,'MaxIter',20)
    addParameter(parser,'rtol',0.05)
    parse(parser,varargin{:})
    
    num_poles = parser.Results.NPoles;
    max_iter = parser.Results.MaxIter;
    rtol = parser.Results.rtol;
      
    if ~isempty(mu_x01)
        mu_resid = data.mu - lm_eval(mu_x01,data.freq);
        [mu_a2,mu_b2] = lm_ord2ratfit(data.freq,mu_resid,num_poles,max_iter,rtol);
    else
        % if empty vector passed for mu_x01, assume mu is fixed at unity.
        % Don't need to estimate poles
        mu_a2 = 0; mu_b2 = 0;
    end
    eps_resid = data.eps - lm_eval(eps_x01,data.freq);
    [eps_a2,eps_b2] = lm_ord2ratfit(data.freq,eps_resid,num_poles,max_iter,rtol);
            
%     mu_a2 = mu_A2(1)*(0.01-1i);
%     mu_b2 = mu_B2(1)*(0.01-1i);
%     disp(mu_B2)
%     % put a1 and b1 in as -1 to avoid errors - these are removed later
%     % anyway
%     mu_x0 = lm_flatten_terms(0,-1,-1,mu_a2,mu_b2);
%     [lb_xmu,ub_xmu] = lm_xbounds(mu_x0,1,1);
%     
%     % remove a1 and b1
%     mu_x0 = [mu_x0(1:2) mu_x0(7:end)];
%     lb_xmu = [lb_xmu(1:2) lb_xmu(7:end)];
%     ub_xmu = [ub_xmu(1:2) ub_xmu(7:end)];
%     disp(mu_x0)
    
    % lsqcurvefit
%     typical_vals = mu_x0;
%     typical_vals(typical_vals==0)=1;
%     options = optimoptions('lsqcurvefit','TypicalX',typical_vals);
%     
%     function F = lmpred(x,xdata)
%         a0 = x(1) + x(2)*1i;
%         a2 = x(3) + x(4)*1i;
%         b2 = x(5) + x(6)*1i;
%         x_in = lm_flatten_terms(a0,0,0,a2,b2);
%         [x_in,np1,np2] = lm_repair_x(x_in,1,1);
%         freq = xdata(1:length(xdata)/2);
%         y_pred = lm_eval(x_in,freq,np1,np2);
%         F = [real(y_pred); imag(y_pred)];
%     end
%         
%     fun = @(x,xdata)lmpred(x,xdata);
%     xdata = repmat(data.freq,2,1);
%     ydata = [real(data.eps); imag(data.eps)];
%     [mu_x_opt,mu_rss] = lsqcurvefit(fun,mu_x0,xdata,ydata,[],[],options);
%     F_opt = lmpred(mu_x_opt,xdata);
    

%     function F = objfunc(x,freq,y_act)
%         a0 = x(1) + x(2)*1i;
%         a2 = x(3) + x(4)*1i;
%         b2 = x(5) + x(6)*1i;
%         x_in = lm_flatten_terms(a0,0,0,a2,b2);
%         [x_in,np1,np2] = lm_repair_x(x_in,1,1);
%         y_pred = lm_eval(x_in,freq,np1,np2);
%         F = sum(abs(y_pred-y_act).^2);
%     end

%     % patternsearch
%     fun = @(x)objfunc(x,data.freq,data.eps);
%     options = optimoptions('patternsearch','MaxIterations',3000,'ScaleMesh',true);
%     [mu_x_opt,mu_rss] = patternsearch(fun,mu_x0,[],[],[],[],lb_xmu,ub_xmu);

%     % particle swarm
%     rng default
%     fun = @(x)objfunc(x,data.freq,data.eps);
%     disp(lb_xmu)
%     disp(ub_xmu)
%     [mu_x_opt,mu_rss,exitflag,output] = particleswarm(fun,length(mu_x0),lb_xmu,ub_xmu);

%     mu_a0 = mu_x_opt(1) + mu_x_opt(2)*1i;
%     mu_a2 = mu_x_opt(3) + mu_x_opt(4)*1i;
%     mu_b2 = mu_x_opt(5) + mu_x_opt(6)*1i;
    
end
    
    