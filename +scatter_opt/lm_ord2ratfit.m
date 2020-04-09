function [A2,B2] = lm_ord2ratfit(freq,actual,varargin)
% Use rationalfit to estimate 2nd-order Laurent model terms. The
% second-order model is given by
%   y = A/(jw-B)^2,
% while rationalfit fits a model of the form
%   F_rf = A/(jw-B).
% Thus, to estimate the A and B terms for the second-order model, we can
% fit rationalfit to the function 
%   F_target = y*(jw-B_est),
% where B_est is our estimate of the appropriate B term. The rationalfit 
% output provides an updated estimate of B_est. By iterating and updating 
% B_est at each iteration, we can converge to a reasonable estimate for B.
% Only one (real) B term can be estimated with this method; additional
% terms can be estimated by estimating the first B and then fitting to the
% residual.
% Parameters
% ----------
% freq : frequencies
% actual: quantity to model
% NPoles : number of 2nd-order poles
% MaxIter: maximum number of iterations to reach convergence
% rtol: relative tolerance for convergence
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'NPoles',1)
    addOptional(parser,'MaxIter',20)
    addOptional(parser,'rtol',0.05)
    parse(parser,varargin{:})
    
    num_poles = parser.Results.NPoles;
    max_iter = parser.Results.MaxIter;
    rtol = parser.Results.rtol;
    
    nfits = rfmodel.rational.empty(num_poles,0);
    A2 = zeros(1,num_poles);
    B2 = zeros(1,num_poles);
    
    % suppress warnings from rationalfit
    warning('off','rf:rationalfit:ErrorToleranceNotMet')
    warning('off','rf:rationalfit:CheckYourData')
    
    for np = 1:num_poles
%         disp(['NPoles=',num2str(np)])
        if np==1
            % First pole fits the function
            resid = actual;
        else
            % Subsequent poles fit the remaining residual
            x_prev = lm_flatten_terms(0,0,0,A2(1:np-1),B2(1:np-1));
            pred = lm_eval(x_prev,freq,1,np-1);
            resid = actual - pred;
        end
        i = 1;
        % initial guess for rationalfit A (pole - B in Laurent model)
        A = -1e8;
        fits = rfmodel.rational.empty(0,0);
        deltas = [];
        while i<=max_iter
            rootdeno= (1j*2*pi*freq - A);
            yfunc = rootdeno.*resid;
            fits(i) = rationalfit(freq,yfunc,'NPoles',1);
            A_new = -sqrt(fits(i).A*A);
%             disp(['Iteration ',num2str(i),': A=',num2str(A,'%10.2e'),' A_new=',...
%                 num2str(A_new,'%10.2e')])
            deltas(i) = abs((A_new - A)/A);
            A = A_new;
            if deltas(i) <= rtol
                converged = true;
                break
            elseif i==max_iter
                converged = false;
            end 
            i = i + 1;
        end
        
        if converged
            % Last fit met rtol
            nfits(np) = fits(end);
        else
            % Last fit did not meet rtol. Find iteration with smallest delta
            [dmin,imin] = min(deltas);
            nfits(np) = fits(imin);
        end
        
        A2(np) = nfits(np).C;
        B2(np) = nfits(np).A;
    end
    
    % turn warnings back on
    warning('on','rf:rationalfit:ErrorToleranceNotMet')
    warning('on','rf:rationalfit:CheckYourData')
end