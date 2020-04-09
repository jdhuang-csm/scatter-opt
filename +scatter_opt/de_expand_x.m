function [p_vec,z_vec,C] = de_expand_x(x)
    n_poles = (length(x)-2)/3;
    p_vec = x(1:n_poles);
    z_vec = x(1+n_poles:2*n_poles) + x(1+2*n_poles:3*n_poles)*1i;
    C = x(end-1) + 1i*x(end);
end