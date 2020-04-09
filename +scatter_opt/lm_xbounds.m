function [lb,ub] = lm_xbounds(x0,num_poles1,num_poles2)
% Get lower and upper bounds for parameter vector. Used in lm_sqfit
% Parameters
% ----------
% x0 : initial parameter vector that has been unpaired

% For real poles (denominator terms): real must be < 0, imag must be 0
% For complex poles: real must be < 0, imag must retain sign
    % specify bounds as complex numbers so that lm_flatten_terms can be
    % used to format them correctly
% For real residues (numerator terms): real must maintain sign, imag must be 0
% For complex residues: real and imag must both maintain sign
    import scatter_opt.*
	[a0,a1,b1,a2,b2] = lm_expand_x(x0,num_poles1,num_poles2);
    
    % no bounds on a0
    a0_lb = -inf - 1i*inf;
    a0_ub = inf + 1i*inf;
    
    % function to get bounds for term vector
    function [alb,aub] = bounds(a_vec,is_opt,is_pole)
    % is_opt: is optimized (i.e. by rationalfit)
    % is_pole: is pole (numerator) vector
        alb = zeros(1,length(a_vec));
        aub = zeros(1,length(a_vec));
        for i=1:length(a_vec)
            ai = a_vec(i);
            % get real bounds first
            if real(ai) > 0
                if is_pole
                    error('Pole lies in right half plane')
                elseif is_opt
                    % if optimized residue, must retain sign
                    alb(i) = 0;
                    aub(i) = inf;
                else
                    % if unoptimized residue, leave unbounded
                    alb(i) = -inf;
                    aub(i) = inf;
                end
            elseif real(ai) < 0
                if is_opt
                    % if optimized, must retain sign
                    alb(i) = -inf;
                    aub(i) = 0;
                else
                    % if unoptimized, leave unbounded
                    alb(i) = -inf;
                    aub(i) = inf;
                end
            else
                % If real part is zero, leave unconstrained
                % This should apply only to 2nd-order residues, which 
                % are initialized with zero real component
                if is_opt
                    error('Optimized term has zero real component')
                elseif is_pole
                    error('Pole has zero real component')
                else
                    alb(i) = -inf;
                    aub(i) = inf;
                end
            end
            % add imag bounds
            if imag(ai) > 0
                if is_opt || is_pole
                    % if optimized or pole, must retain sign
                    alb(i) = alb(i) + 1i*1e-5;
                    aub(i) = aub(i) + 1i*inf;
                else
                    % if unoptimized residue, leave unconstrained
                    alb(i) = alb(i) - 1i*inf;
                    aub(i) = aub(i) + 1i*inf;
                end
            elseif imag(ai) < 0
                if is_opt
                    % if optimized, must retain sign
                    alb(i) = alb(i) - 1i*inf;
                    aub(i) = aub(i) + 1i*0;
                elseif is_pole
                    error('Pole has negative imaginary component')
                else
                    % if unoptimized residue, leave unconstrained
                    alb(i) = alb(i) - 1i*inf;
                    aub(i) = aub(i) + 1i*inf;
                end
            else
                % if real, imag==0 (i.e. must stay real)
                alb(i) = alb(i) + 1i*0;
                aub(i) = aub(i) + 1i*0;
            end
        end
    end
    
    % get bounds for all poles
    [a1_lb,a1_ub] = bounds(a1,true,false);
    [b1_lb,b1_ub] = bounds(b1,true,true);
    [a2_lb,a2_ub] = bounds(a2,false,false);
    [b2_lb,b2_ub] = bounds(b2,false,true);
    
    % flatten complex arrays
    lb = lm_flatten_terms(a0_lb,a1_lb,b1_lb,a2_lb,b2_lb);
    ub = lm_flatten_terms(a0_ub,a1_ub,b1_ub,a2_ub,b2_ub);
end