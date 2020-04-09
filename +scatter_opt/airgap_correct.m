function [mu_c,eps_c] = airgap_correct(mu,eps,b,gap)
    d = b-gap;
    eps_c_re = real(eps)*d./(b-(b-d)*real(eps));
    eps_c_im = eps_c_re.*(imag(eps)./real(eps)).*(b./(b-(b-d).*real(eps)));
    eps_c = eps_c_re + 1i*eps_c_im;
    
    mu_c_re = real(mu)*(b/d)-(b-d)/d;
    mu_c_im = imag(mu)*b/d;
    mu_c = mu_c_re + 1i*mu_c_im;
end
    