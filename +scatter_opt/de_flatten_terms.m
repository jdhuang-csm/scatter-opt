function x = de_flatten_terms(p_vec,z_vec,C)
    n_poles = length(p_vec);    
    x = zeros(1,3*n_poles+2);
    x(1:n_poles) = p_vec;
    x(1+n_poles:2*n_poles) = real(z_vec);
    x(1+2*n_poles:3*n_poles) = imag(z_vec);
    x(end-1) = real(C);
    x(end) = imag(C);
end