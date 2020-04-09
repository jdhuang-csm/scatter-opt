function axes = plot_PP(data,varargin)
% plot mu and epsilon
    import scatter_opt.*
	parser = inputParser;
    addOptional(parser,'axes',{})
    addParameter(parser,'VariableNames',{'mu' 'eps'})
    addParameter(parser,'label','')
    parse(parser,varargin{:})
    
    if isempty(parser.Results.axes)
        ax1 = subplot(1,2,1); hold on
        ax2 = subplot(1,2,2); hold on
        axes = {ax1 ax2};
    else
        axes = parser.Results.axes;
        ax1 = axes{1}; hold(ax1,'on')
        ax2 = axes{2}; hold(ax2,'on')
    end
    
    label = parser.Results.label;

    mu_name = parser.Results.VariableNames{1};
    eps_name = parser.Results.VariableNames{2};
    plot(ax1,data.freq,real(data.(mu_name)),'DisplayName',strcat(label,' Re'))
    plot(ax1,data.freq,-imag(data.(mu_name)),'DisplayName',strcat(label,' Im'))
    title(ax1,'Mu')
    legend(ax1)
    hold(ax1,'off')

    plot(ax2,data.freq,real(data.(eps_name)),'DisplayName',strcat(label,' Re'))
    plot(ax2,data.freq,-imag(data.(eps_name)),'DisplayName',strcat(label,' Im'))
    title(ax2,'Epsilon')
    hold(ax2,'off')
    
    for i=1:2
        ax = axes{i};
        grid(ax,'on')
        xlabel(ax,'Frequency (Hz)')
        legend(ax,'Orientation','horizontal','Location','southoutside')
    end
    
end