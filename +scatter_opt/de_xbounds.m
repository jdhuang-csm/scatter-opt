function [lb,ub] = de_xbounds(x)
    n_poles = (length(x)-2)/3;
    lb = zeros(1,length(x));
    ub = zeros(1,length(x));
    % poles must be negative
    lb(1:n_poles) = -inf;
    ub(1:n_poles) = 0;
    % zeros and C are unconstrained
    lb(1+n_poles:end) = -inf;
    ub(1+n_poles:end) = inf;
end
    