function y = lm_eval(x,freq,varargin)
% Evalute 2nd-order Laurent model
% Paremeters
% ----------
% x : parameter vector
% freq : frequencies to evaluate
% num_poles1 : number of 1st-order poles. If not specified, assumes that
%   num_poles1=num_poles2=(length(x)-2)/8
% num_poles1 : number of 2nd-order poles
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'num_poles1',-1)
    addOptional(parser,'num_poles2',-1)
    parse(parser,varargin{:})
    if isempty(x)
        y = ones(length(freq),1);
    else
        if parser.Results.num_poles1==-1
            % if num_poles not specified, assume same number of 1st- and
            % 2nd-order poles
            num_poles1 = round((length(x)-2)/8);
            num_poles2 = num_poles1;
            if 2+(num_poles1+num_poles2)*4~=length(x)
                error('Different number of 1st- and 2nd-order poles! Specify num_poles1 and num_poles2.')
            end
        elseif parser.Results.num_poles1~=-1 && parser.Results.num_poles2~=-1
            % Number of 1st- and 2nd-order poles specified in function call
            num_poles1 = parser.Results.num_poles1;
            num_poles2 = parser.Results.num_poles2;
        else
            error('Both num_poles1 and num_poles2 must be specified')
        end


        [a0,a1,b1,a2,b2] = lm_expand_x(x,num_poles1,num_poles2);
        omega = 2.*pi.*freq;
        % Matrix of denominators for 1st-order terms
        % freq is Nx1, b1 is 1xp, freq*b1 is Nxp
        deno_mat1 = 1i*omega-b1;
        frac_mat1 = a1./deno_mat1;
        sum1 = sum(frac_mat1,2);

        % Matrix of denominators for 2nd-order terms
        deno_mat2 = (1i*omega-b2).^2; %(1+1i*omega*b2).^2;
        frac_mat2 = a2./deno_mat2;
        sum2 = sum(frac_mat2,2);

        y = a0 + sum1 + sum2;
    end
end
    
    
    