function lm_plotfit(x,actual,freq,varargin)
% Plot complex Laurent model vs. frequency
% Paremeters
% ----------
% x : parameter vector
% freq : frequencies to evaluate
% num_poles1 : number of 1st-order poles. If not specified, assumes that
% num_poles1=num_poles2=(length(x)-2)/8
% num_poles1 : number of 2nd-order poles
    import scatter_opt.*
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
    
    axr = subplot(1,2,1);
    axi = subplot(1,2,2);
    pred = lm_eval(x,freq,num_poles1,num_poles2);
    plot(axr,freq,real(actual),freq,real(pred))
    plot(axi,freq,imag(actual),freq,imag(pred))
    legend(axr,'Actual','Predicted')
    xlabel(axr,'Frequency (Hz)')
    ylabel(axr,'Real')
    xlabel(axi,'Frequency (Hz)')
    ylabel(axi,'Imag')
end