function points = lm_mspoints(x0,num_points,shift_factor)
% Generate random starting points for MultiStart
% Parameters
% ----------
% x0 : initial parameter vector
% num_points : number of points to generate
% shift_factor: factor by which parameters may be multiplied or divided to
%   generate random points
    rng default
    nvar = length(x0);
    % uniformly spaced values in [0,1)
    rands = rand(num_points,nvar);
    % transform linear-uniformly spaced values in [0,1) range to 
    % logarithmic-uniformly spaced values in [1/shift_factor,shift_factor) 
    % range
    factors = (1/shift_factor)*exp(rands*2*log(shift_factor));
    % multiply original coordinates by factors to get shifted coordinates
    points = factors.*x0;
end