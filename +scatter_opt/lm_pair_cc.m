function [a_cc,b_cc] = lm_pair_cc(a_vec,b_vec)
% Identify complex poles in Laurent model term vectors and add their complex
% conjugates to the vectors
% Parameters
% ----------
% a_vec : A-vector of residue terms for Laurent model (i.e. a1 or a2)
% b_vec : B-vector of pole terms for Laurent model (i.e. b1 or b2)
% add_existing : if true, add complex conjugates to vectors even if they
%   appear to already exist. This is necessary when running lsqcurvefit, as
%   the optimizer may shift some terms such that they are nearly 
    % order matters! need to keep order of a_vec and b_vec aligned
    
    if length(a_vec)~=length(b_vec)
        error('a_vec and b_vec must be same length')
    end
    a_cc = a_vec;
    b_cc = b_vec;
    for i=1:length(a_vec)
        ai = a_vec(i);
        bi = b_vec(i);
        if imag(bi) > 0 
            if ~ismember(round(conj(ai),10),round(a_vec,10)) &&...
                    ~ismember(round(conj(bi),10),round(b_vec,10))
                % If ai has positive imag part, add its complex conjugate
                % and the corresponding conjugate for bi
                if imag(ai)~=0
                    a_cc = [a_cc conj(ai)];
                    b_cc = [b_cc conj(bi)];
                else
                    err = ['bi is complex but ai is real. ai=',...
                        num2str(ai),', bi=',num2str(bi)];
                    error(err)
                end
            elseif imag(bi) < 1e-5 && imag(ai)==0
                % This is actually supposed to be a real pole, but the
                % solver shifted it by a small amount away from zero. No
                % complex conjugate should be added
                continue
            else
                % If the complex conjugate to ai or bi is already in a_vec
                % or b_vec, throw an error
                disp(i)
                disp(ai)
                disp(round(conj(ai),10)==round(a_vec,10))
                disp(bi)
                disp(round(conj(bi),10)==round(b_vec,10))
                error('Complex conjugate already exists in a_vec or b_vec')
            end
        elseif imag(bi) < 0
            % there shouldn't be any negative imaginary components in a_vec
            disp(bi)
            error('Unexpected negative imaginary component in bi')
        else
            % bi is real. Don't need to add complex conjugate
            if imag(ai) < 1e-5 % allow some margin for solver
                continue
            else
                err = ['bi is real but ai is complex. ai=',...
                    num2str(ai),', bi=',num2str(bi)];
                error(err)
                
            end
        end
    end
%     cc = conj(a_vec);
%     a_cc = [a_vec cc];
%     a_cc = unique(a_cc);
end
            
            
