function N = nrw_refracIndex(Sdata,L,varargin)
% Calculate complex refractive index for KK determination of branch
% Parameters
% ----------
% Sdata: table with freq, s11, and s21
% L : sample thickness
% Named Parameters
% ----------------
% Branch: which branch to use. Affects real part only
	import scatter_opt.*
    parser = inputParser;
    addOptional(parser,'Branch',0)
    parse(parser,varargin{:})
    
    p = parser.Results.Branch;
    omega = Sdata.freq*2*pi;

    % impedance
    Z = nrw_eta(Sdata.s11,Sdata.s21,'Relative',false); 
    % Re(Z) >= 0
    Zr_sign = sign(real(Z));
    % if real part is exactly zero, leave sign
    Zr_sign(Zr_sign==0) = 1;
    Z = Z.*Zr_sign;
    
    R01 = (Z - 1)./(Z + 1);
    ei = Sdata.s21./(1 - Sdata.s11.*R01);
    
    % free space wavenumber
    c_vac = 299792458.000176;
    k0 = omega./c_vac;
    
    kap = -real(log(ei))./(k0*L);
    n = imag(log(ei))./(k0*L) + 2*p*pi./(k0*L);
    N = n + 1i*kap;
end
    
    