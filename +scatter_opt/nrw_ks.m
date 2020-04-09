function ks = nrw_ks(s11,s21,L,p)
% Wavenumber from Arslanagic et al. 2013 (eq. 24)
% Parameters:
% -----------
% s11: complex s11 values
% s21: complex s21 values
% L : sample thickness (d in paper)
% p : branch of argument of Z (int)
	import scatter_opt.*
    Z = nrw_Z(s11,s21);
    ks = (1./L).*(-(angle(Z)+2.*pi.*p) + 1i.*log(abs(Z)));
end