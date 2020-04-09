function nkk = nrw_refracKK(freq,kap)
% calculate real part of refractive index from imag part using
% Kramers-Kronig relations. Method from:
% Z. Szabó, G. H. Park, R. Hedge, and E. P. Li, “A unique extraction of metamaterial parameters based on Kramers-Kronig relationship,” IEEE Trans. Microw. Theory Tech., vol. 58, no. 10, pp. 2646–2653, 2010, doi: 10.1109/TMTT.2010.2065310.
% Parameters
% ----------
% freq: frequency values
% kap: imaginary part of refractive index
    omega = 2*pi*freq;
    
    function f = intfun(x,omega,kap)
        % integrand function
        f = (omega.*kap)./(omega.^2 - x.^2);
    end

    function int = cauchy(x,omega,kap)
        idx = find(omega==x,1);
        fint = intfun(x,omega,kap);
        if idx<=2
            int = trapz(omega(idx:end),fint(idx:end));
        elseif idx>=length(omega)-1
            int = trapz(omega(1:idx-1),fint(1:idx-1));
        else
            int = trapz(omega(1:idx-1),fint(1:idx-1))...
                + trapz(omega(idx+1:end),fint(idx+1:end));
        end
    end

    fun = @(x)cauchy(x,omega,kap);
    nkk = (1 + (2/pi)*arrayfun(fun,omega));
%     bkk = fun(omega(1));
end