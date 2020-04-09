function [x,num_poles1,num_poles2] = lm_prep_x0(x0,num_ccp,num_poles1,num_poles2)
% Prepare initial parameter vector for lsq fit. Unpair a1 and b1.
% Initialize specified number of cc pairs in a2 and b2, and then unpair.
% Parameters
% ----------
% x0 : initial guess parameter vector with a0, a1, and b1 extracted from
%   rationalfit. Must contain CC pairs (can't be unpaired)
% num_ccp : number of complex conjugate pairs to initialize in a2 and b2
% EstOrd2 : if true, estimate 2nd-order pole starting values
    import scatter_opt.*
	[a0,a1,b1,a2,b2] = lm_expand_x(x0,num_poles1,num_poles2);
    num_poles = length(a2);
    if num_ccp>floor(num_poles/2)
        error('Too many complex conjugate pairs')
    end
    % unpair complex conjugates in a1 and b1
    [a1,b1] = lm_unpair_cc(a1,b1);
    
    if sum(a2)==0 && sum(b2)==0
        % If a2 and b2 are not set, initialize small imag components for 
        % poles that will be part of complex conjugate pairs
        for i=1:num_ccp
            % for each cc pair, initialize pair for a2 and b2
            % numerator (residue): initialize small imag component to indicate
            % complex. Sign will not necessarily be maintained during
            % optimization
            a2(i) = 1e10*(1 + 1i);
            a2(i+num_ccp) = conj(a2(i));
            % denominator (pole) must lie in left half-plane. Imag sign will be
            % fixed during optimization, since the complex conjugate has the
            % opposite sign
            b2(i) = 1e8*(-1 + 1i);
            b2(i+num_ccp) = conj(b2(i));
        end

        for i=2*num_ccp+1:num_poles
            % initialize real 2nd-order residues at 0 (code explicitly for
            % clarity)
            a2(i) = 1e10;
            % Initialize real 2nd-order poles as small negative numbers
            b2(i) = -1e8;
        end
    else
        for i=1:num_ccp
            % assuming a2 and b2 were estimated as negative real numbers, 
            % convert to imaginary
            a2(i) = a2(i)*(0.1 - 1i);
            a2(i+num_ccp) = conj(a2(i));
            % denominator (pole) must lie in left half-plane. Imag sign will be
            % fixed during optimization, since the complex conjugate has the
            % opposite sign
            b2(i) = b2(i)*(0.1 + 1i);
            b2(i+num_ccp) = conj(b2(i));
        end
        
    end
        
    [a2,b2] = lm_unpair_cc(a2,b2);
    num_poles1 = length(a1);
    num_poles2 = length(a2);
    if num_poles1~=length(b1)
        err = ['Different number of terms in A1 and B1.\nA1: ',...
            num2str(a1),'\nB1: ',num2str(b1)];
        error(err)
    elseif num_poles2~=length(b2)
        err = ['Different number of terms in A2 and B2.\nA2: ',...
            num2str(a2),'\nB2: ',num2str(b2)];
        error(err)
    end
    x = lm_flatten_terms(a0,a1,b1,a2,b2);
end