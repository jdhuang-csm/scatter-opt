function out = de_output(defit,Sdata,PPdata,varargin)
% Output full result set with option to save to file
% Parameters
% ----------
% lmfit : lmfit object obtained frmo lm_lsqfit
% Sdata : data table with measured freq and s_ij
% PPdata: data table with freq and extracted mu and eps
% Named (optional) parameters
% ---------------------------
% SaveFile : file to save output to
% ExtractMethod : method of extraction ('nrw' or 'nni'). If both 
%   ExtractMethod and SaveFile are specified, add ExtractMethod to end of 
%   mu and eps column names in output file
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'SaveFile','')
    addOptional(parser,'ExtractMethod','')
    parse(parser,varargin{:})
    
    mu_pred = ones(height(Sdata),1);
    eps_pred = de_eval(defit.x_eps,Sdata.freq);
    s_pred = de_rev_eval(defit,Sdata.freq);
    
    out = Sdata;
    
    save_file = parser.Results.SaveFile;
    extract_method = parser.Results.ExtractMethod;
    if ~strcmp(save_file,'') && ~strcmp(extract_method,'')
        % if SaveFile and ExtractMethod given, add suffix to indicate
        % extraction method
        out.(strcat('mu_',extract_method)) = PPdata.mu;
        out.(strcat('eps_',extract_method)) = PPdata.eps;
    else
        out.mu = PPdata.mu;
        out.eps = PPdata.eps;
    end
    out.s11_model = s_pred.s11;
    out.s21_model = s_pred.s21;
    out.s12_model = s_pred.s12;
    out.s22_model = s_pred.s22;
    out.mu_model = mu_pred;
    out.eps_model = eps_pred;
    
    if ~strcmp(save_file,'')
        writetable(out,save_file)
    end
end