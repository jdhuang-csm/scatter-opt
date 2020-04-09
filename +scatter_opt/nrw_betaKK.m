function bkk = nrw_betaKK(freq,alpha)
% Estimate imag part of propagation constant from real part using KK relations
	import scatter_opt.*
    omega = 2*pi*freq;
    c_vac = 299792458.000176;
    % free space wavenumber
    b0 = omega./c_vac;
    
    function f = intfun(x,omega,alpha,b0)
        % integrand function
        f = (omega.*alpha./b0)./(omega.^2 - x.^2);
    end

    function int = cauchy(x,omega,alpha,b0)
        idx = find(omega==x,1);
        fint = intfun(x,omega,alpha,b0);
        if idx<=2
            int = trapz(omega(idx:end),fint(idx:end));
        elseif idx>=length(omega)-1
            int = trapz(omega(1:idx-1),fint(1:idx-1));
        else
            int = trapz(omega(1:idx-1),fint(1:idx-1))...
                + trapz(omega(idx+1:end),fint(idx+1:end));
        end
    end

    fun = @(x)cauchy(x,omega,alpha,b0);
    bkk = b0.*(1 + (2/pi)*arrayfun(fun,omega));
%     bkk = (2/pi)*arrayfun(fun,omega);
end
    