function gam = rev_gamma(omega,mu,eps,lambda_c)
    % propagation constant in material
    c_vac = 299792458.000176;
    gam = 1i.*sqrt((omega.^2.*mu.*eps)./(c_vac.^2) - (2.*pi./lambda_c).^2);
end