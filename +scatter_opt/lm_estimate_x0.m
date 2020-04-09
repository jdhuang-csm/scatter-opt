function [mu_x0,eps_x0,mu_np1,eps_np1,np2] = lm_estimate_x0(nrw_out,varargin)
% Estimate initial parameter vectors for mu and eps including first- and
% second-order poles.
% Parameters
% ----------
% nrw_out : table with freq and extracted mu and epsilon
% NPoles1: number of 1st-order poles to use for both mu and eps
% NPoles2: number of 2nd-order poles to use for both mu and eps
% EstOrd2: if true, estimate initial parameters for 2nd-order poles. If
%   false, initialize 2nd-order poles parameters at zero
% MuNPoles1: number of 1st-order poles to use for mu. If specified,
%   overrides NPoles1
% EpsNPoles1: number of 1st-order poles to use for eps. If specified,
%   overrides NPoles1
    import scatter_opt.*
	parser = inputParser;
    addParameter(parser,'NPoles1',2)
    addParameter(parser,'NPoles2',1)
    addParameter(parser,'EstOrd2',true)
    addParameter(parser,'MuNPoles1',-1)
    addParameter(parser,'EpsNPoles1',-1)
    parse(parser,varargin{:})
    
    % get 1st-order fits
    if parser.Results.MuNPoles1~=-1
        mu_np1 = parser.Results.MuNPoles1;
    else
        mu_np1 = parser.Results.NPoles1;
    end
    if parser.Results.EpsNPoles1~=-1
        eps_np1 = parser.Results.EpsNPoles1;
    else
        eps_np1 = parser.Results.NPoles1;
    end
    if mu_np1==eps_np1
        init_result = lm_initialfit(nrw_out,'NPoles',mu_np1);
        mu_fit = init_result.mu_fit(1);
        eps_fit = init_result.eps_fit(1);
    elseif mu_np1==0
        eps_result = lm_initialfit(nrw_out,'NPoles',eps_np1);
        eps_fit = eps_result.eps_fit(1);
    else
        mu_result = lm_initialfit(nrw_out,'NPoles',mu_np1);
        mu_fit = mu_result.mu_fit(1);
        eps_result = lm_initialfit(nrw_out,'NPoles',eps_np1);
        eps_fit = eps_result.eps_fit(1);
    end
    
    if mu_np1==0
       mu_x01 = [];
    else
        mu_x01 = lm_extract_rfx(mu_fit);
        [mu_a0,mu_a1,mu_b1,mu_a2,mu_b2] = lm_expand_x(mu_x01);
    end
    
    eps_x01 = lm_extract_rfx(eps_fit);
    [eps_a0,eps_a1,eps_b1,eps_a2,eps_b2] = lm_expand_x(eps_x01);
    
    if parser.Results.EstOrd2
        % estimate 2nd-order terms
%         [mu_a2,mu_b2] = lm_ord2ratfit(nrw_out.freq,nrw_out.mu);
%         [eps_a2,eps_b2] = lm_ord2ratfit(nrw_out.freq,nrw_out.eps);
        [mu_a2,mu_b2,eps_a2,eps_b2] = lm_initialfit2(nrw_out,mu_x01,...
            eps_x01,'NPoles',parser.Results.NPoles2);
%         if parser.Results.NPoles2 > 1
%             mu_a2 = repmat(mu_a2,1,parser.Results.NPoles2);
%             mu_b2 = repmat(mu_b2,1,parser.Results.NPoles2);
%             eps_a2 = repmat(eps_a2,1,parser.Results.NPoles2);
%             eps_b2 = repmat(eps_b2,1,parser.Results.NPoles2);
%         end
    else
        % set 2nd-order terms to zero
        mu_a2 = zeros(parser.Results.NPoles2,1);
        mu_b2 = zeros(parser.Results.NPoles2,1);
        eps_a2 = zeros(parser.Results.NPoles2,1);
        eps_b2 = zeros(parser.Results.NPoles2,1);
    end
    
    np2 = parser.Results.NPoles2;
    
    if mu_np1==0
        mu_x0 = [];
    else
        mu_x0 = lm_flatten_terms(mu_a0,mu_a1,mu_b1,mu_a2,mu_b2);
    end
    eps_x0 = lm_flatten_terms(eps_a0,eps_a1,eps_b1,eps_a2,eps_b2);
end
        
    