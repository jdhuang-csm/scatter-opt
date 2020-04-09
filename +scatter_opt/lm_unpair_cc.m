function [a_up,b_up] = lm_unpair_cc(a_vec,b_vec)
% Identify complex conjugate pairs in Laurent model term vectors and remove
% the poles with negative imag components. Also remove the corresponding
% residues, which may have different signs of the imag component.
% Parameters
% ----------
% a_vec : A-vector of residues for Laurent model (i.e. A1 or A2)
% b_vec : B-vector of poles for Laurent model (i.e. B1 or B2)
    if length(a_vec)~=length(b_vec)
        error('a_vec and b_vec must be same length')
    end
    a_up = [];
    b_up = [];
    for i=1:length(a_vec)
        ai = a_vec(i);
        bi = b_vec(i);
        if imag(bi) < 0 
            if ismember(round(conj(ai),10),round(a_vec,10)) &&...
                    ismember(round(conj(bi),10),round(b_vec,10))
                % If ai is part of a complex conjugate pair, don't add ai 
                % or the corresponding bi to the output vector 
                continue
            else
                % If the complex comjugate to ai or bi is missing, raise a
                % warning and add conj(ai) and conj(bi) to the output
                % vectors
                warning('Complex conjugate missing in a_vec or b_vec')
                a_up = [a_up conj(ai)];
                b_up = [b_up conj(bi)];
            end
        else
            % if bi is a real number or a complex number with positive 
            % imaginary part, put bi and the corresponding ai in the output
            % vectors
            a_up = [a_up ai];
            b_up = [b_up bi];
        end
    end
end