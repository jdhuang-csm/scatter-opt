function Z = nrw_Z(s11,s21)
% Z parameter (method A) from Arslanagic et al. 2013 (eq. 23a)
    mu_0 = 1.25663706e-6;
    eps_0 = 8.85418782e-12;
    % intrinsic impedance of free space
    eta_0 = sqrt(mu_0/eps_0);
    % get impedance
    eta = nrw_eta(s11,s21);
    
    nume = s21.*(eta+eta_0);
    deno = (eta+eta_0) - s11.*(eta-eta_0);
    Z = nume./deno;
end
    
    