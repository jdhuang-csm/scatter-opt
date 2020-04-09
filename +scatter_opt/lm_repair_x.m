function [x_out,np1_out,np2_out] = lm_repair_x(x,num_poles1,num_poles2)
% Re-pair unpaired complex conjugates 
% Parameters
% ---------
% x : unpaired parameter vector
% num_poles1 : number of 1st-order poles in x
% num_poles2 : number of 2nd-order poles in x
    import scatter_opt.*
	[a0,a1,b1,a2,b2] = lm_expand_x(x,num_poles1,num_poles2);
    [a1,b1] = lm_pair_cc(a1,b1);
    [a2,b2] = lm_pair_cc(a2,b2);
    np1_out = length(a1);
    np2_out = length(a2);
    if np1_out~=length(b1)
        err = ['Different number of terms in A1 and B1.\nA1: ',...
            num2str(a1),'\nB1: ',num2str(b1)];
        error(err)
    elseif np2_out~=length(b2)
        err = ['Different number of terms in A2 and B2.\nA2: ',...
            num2str(a2),'\nB2: ',num2str(b2)];
        error(err)
    end
    x_out = lm_flatten_terms(a0,a1,b1,a2,b2);
end