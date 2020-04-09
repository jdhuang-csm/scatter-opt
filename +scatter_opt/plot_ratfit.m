function plot_ratfit(fit,actual,freq,varargin)
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'FlipIm',false)
    addOptional(parser,'axes',{})
    parse(parser,varargin{:})
    
    if isempty(parser.Results.axes)
        axr = subplot(1,2,1);
        axi = subplot(1,2,2);
    else
        axr = parser.Results.axes{1}; hold(axr,'on')
        axi = parser.Results.axes{2}; hold(axi,'on')
    end
    
    [pred,outfreq] = freqresp(fit,freq);
    
    if parser.Results.FlipIm
        % flip sign of imag part to meet sign convention in plot
        actual = conj(actual);
        pred = conj(pred);
    end
    
    plot(axr,freq,real(actual),freq,real(pred))
    plot(axi,freq,imag(actual),freq,imag(pred))
    legend(axr,'Actual','Predicted')
    xlabel(axr,'Frequency (Hz)')
    ylabel(axr,'Real')
    xlabel(axi,'Frequency (Hz)')
    ylabel(axi,'Imag')
    
    grid(axr,'on')
    grid(axi,'on')
    
    hold(axr,'off')
    hold(axi,'off')
end