function axes = plot_rivf(freq,z,varargin)
    % args: freq, values, real axis, imag axis
    % plot real and imaginary parts vs. frequency
    parser = inputParser;
    addOptional(parser,'axes',{})
    addParameter(parser,'label','')
    parse(parser,varargin{:})
    
    if isequal(parser.Results.axes,{})
        % no axes provided
        axr = subplot(1,2,1);
        axi = subplot(1,2,2);
        axes = {axr,axi};
    else
        % existing axes provided
        axes = parser.Results.axes;
        axr = axes{1};
        axi = axes{2};
        % turn hold on to not overwrite existing lines
        hold(axr,'on')
        hold(axi,'on')
    end
    %{
    if nargin-2==2
        axr = cell2mat(varargin(1));
        disp(axr)
        axi = varargin(2);
        axes 
    end
    %}
    plot(axr,freq,real(z),'DisplayName',parser.Results.label)
    plot(axi,freq,imag(z),'DisplayName',parser.Results.label)
    xlabel(axr,'Frequency (Hz)')
    ylabel(axr,'Real')
    xlabel(axi,'Frequency (Hz)')
    ylabel(axi,'Imag')
    
    % turn hold off
    hold(axr,'off')
    hold(axi,'off')
    
    grid(axr,'on')
    grid(axi,'on')
end