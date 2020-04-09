function eta = nrw_eta(s11,s21,varargin)
% impedance from Arslanagic et al. 2013 (eq. 22)
    parser = inputParser;
    addOptional(parser,'Relative',true)
    parse(parser,varargin{:})
    
    % impedance of sample
    nume = (s11+1).^2 - s21.^2;
    deno = (s11-1).^2 - s21.^2;
    eta = sqrt(nume./deno);
    
    if parser.Results.Relative
        % normalize for impedance of free space
        mu_0 = 1.25663706e-6;
        eps_0 = 8.85418782e-12;
        % intrinsic impedance of free space
        eta_0 = sqrt(mu_0/eps_0);
        eta = eta_0.*eta;
    end
end