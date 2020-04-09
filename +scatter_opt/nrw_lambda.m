function lam = nrw_lambda(tc,L,varargin)
% Big lambda (Annex 2 eq. 1.5)
% Parameters
% ----------
% tc: transmission coefficient
% L : sample thickness in m
% branch: integer, default 0
	import scatter_opt.*
    parser = inputParser;
    addOptional(parser,'branch',0)
    addOptional(parser,'PhaseUnwrap',true)
    parse(parser,varargin{:})
    p = parser.Results.branch;
    
    % calculate log of 1/T based on branch
    beta = angle(1./tc);
    if parser.Results.PhaseUnwrap
        buw = unwrap(beta);
        step_index = find(buw~=beta,1);
        if ~isempty(step_index)
            warning(['Phase unwrap identified branch change at index ',mat2str(step_index)])
        end
        beta = buw;
    end
    log_tc = log(abs(1./tc)) + 1i.*(beta + 2.*pi.*p);
    
    lam = sqrt(-1./((1./(2.*pi.*L)).*log_tc).^2);
end