function [a0,a1,b1,a2,b2] = lm_expand_x(x,varargin)
% Expand parameter vector x into terms for Laurent model. Complex
% conjugates must be paired if num_poles1 and num_poles2 are not specified,
% as unpaired CCs will throw off indexing.
    parser = inputParser;
    addOptional(parser,'num_poles1',-1)
    addOptional(parser,'num_poles2',-1)
    parse(parser,varargin{:})
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
    end
    a0 = x(1) + 1i*x(2);
    a1 = x(3:3+num_poles1-1) + 1i*x(3+num_poles1:3+2*num_poles1-1);
    b1 = x(3+2*num_poles1:3+3*num_poles1-1)...
        + 1i*x(3+3*num_poles1:3+4*num_poles1-1); 
    a2 = x(3+4*num_poles1:3+4*num_poles1+num_poles2-1) ...
        + 1i*x(3+4*num_poles1+num_poles2:3+4*num_poles1+2*num_poles2-1);
    b2 = x(3+4*num_poles1+2*num_poles2:3+4*num_poles1+3*num_poles2-1)...
        + 1i*x(3+4*num_poles1+3*num_poles2:3+4*num_poles1+4*num_poles2-1);    
end