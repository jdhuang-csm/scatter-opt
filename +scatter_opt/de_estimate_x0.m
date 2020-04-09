function x0 = de_estimate_x0(nrw_out,varargin)
% Estimate starting parameters to fit epsilon.
% Parameters
% ----------
% nrw_out : table with freq and extracted mu and eps
% Named (optional) parameters
% ---------------------------
% NPoles : Number of poles to be distributed uniformly acress frequency 
%   range. Defaults to 2
% PoleFreq: pole frequencies. If specified, NPoles is ignored, and one pole
%   is initialized at each frequency specified.
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'NPoles',2)
    addOptional(parser,'PoleFreq',[])
    parse(parser,varargin{:})
    
    [p_vec,z_vec,C_vec] = de_initialfit(nrw_out,parser.Results.NPoles,...
        parser.Results.PoleFreq);
    C = prod(C_vec);
    x0 = de_flatten_terms(p_vec,z_vec,C);
end