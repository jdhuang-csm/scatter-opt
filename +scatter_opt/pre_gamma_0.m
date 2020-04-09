function gam = pre_gamma_0(freq,lambda_c)
% gamma_0 for preprocessing
% turns out to be equivalent to rev_gamma_0
    omega = freq*2*pi;
    kc = 2*pi/lambda_c;
    c_lab = 299704644.54;
    k0 = omega./c_lab;
    gam = sqrt(kc.^2 - k0.^2);
end
    