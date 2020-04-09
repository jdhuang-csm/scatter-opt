function gam = rev_gamma_0(omega,lambda_c)
    % propagation constant in air
    c_lab = 299704644.54;
    gam = 1i.*sqrt((omega./c_lab).^2-(2.*pi./lambda_c).^2);
end