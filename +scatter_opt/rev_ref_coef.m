function rc = rev_ref_coef(mu,gam0,gam)
    % reflection coefficient
    %mu_0 = 1.25663706e-6;
    ratio = mu.*gam0./gam;
    nume = ratio - 1;
    deno = ratio + 1;
    rc = nume./deno;
end