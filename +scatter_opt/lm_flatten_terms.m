function x = lm_flatten_terms(a0,a1,b1,a2,b2)
% Flatten Laurent model terms into parameter vector. Complex conjugates can
% be paired or unpaired, but the a1 and b1 must have the same number of
% complex poles, and a2 and b2 must have the same number of complex poles.
    % first-order poles
    num_poles1 = length(a1);
    % second-order poles
    num_poles2 = length(a2);
    x = zeros(1,2+4*num_poles1+4*num_poles2);
    x(1) = real(a0);
    x(2) = imag(a0);
    x(3:3+num_poles1-1) = real(a1);
    x(3+num_poles1:3+2*num_poles1-1) = imag(a1);
    x(3+2*num_poles1:3+3*num_poles1-1) = real(b1);
    x(3+3*num_poles1:3+4*num_poles1-1) = imag(b1);
    x(3+4*num_poles1:3+4*num_poles1+num_poles2-1) = real(a2);
    x(3+4*num_poles1+num_poles2:3+4*num_poles1+2*num_poles2-1) = imag(a2);
    x(3+4*num_poles1+2*num_poles2:3+4*num_poles1+3*num_poles2-1) = real(b2);
    x(3+4*num_poles1+3*num_poles2:3+4*num_poles1+4*num_poles2-1) = imag(b2);
end