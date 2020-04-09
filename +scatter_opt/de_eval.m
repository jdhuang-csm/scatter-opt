function pred = de_eval(x,freq)
    import scatter_opt.*
	[p_vec,z_vec,C] = de_expand_x(x);
    n_poles = length(p_vec);
    omega = 2*pi*freq;
    pred = C.*ones(length(freq),1);
    for n=1:n_poles
        frac = (1i.*omega - z_vec(n))./(1i.*omega - p_vec(n));
        pred = pred.*frac;
    end
end
        
        