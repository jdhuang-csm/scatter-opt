function axes = plot_Sij(data,varargin)
% plot all S parameters
    parser = inputParser;
    addOptional(parser,'axes',{})
    addParameter(parser,'label','')
    parse(parser,varargin{:})
    
    if isempty(parser.Results.axes)
        ax11 = subplot(2,2,1);
        hold on
        ax21 = subplot(2,2,2);
        hold on
        ax12 = subplot(2,2,3);
        hold on
        ax22 = subplot(2,2,4);
        hold on
        axes = {ax11 ax21 ax12 ax22};
    else
        axes = parser.Results.axes;
        ax11 = axes{1}; hold(ax11,'on')
        ax21 = axes{2}; hold(ax21,'on')
        ax12 = axes{3}; hold(ax12,'on')
        ax22 = axes{4}; hold(ax22,'on')
    end
    
    label = parser.Results.label;
    
    plot(ax11,data.freq,real(data.s11),'DisplayName',strcat(label,' Re'))
    plot(ax11,data.freq,imag(data.s11),'DisplayName',strcat(label,' Im'))
    title(ax11,'S11')
    legend(ax11)
    hold(ax11,'off')

    plot(ax21,data.freq,real(data.s21),'DisplayName',strcat(label,' Re'))
    plot(ax21,data.freq,imag(data.s21),'DisplayName',strcat(label,' Im'))
    title(ax21,'S21')
    hold(ax21,'off')
    
    plot(ax12,data.freq,real(data.s12),'DisplayName',strcat(label,' Re'))
    plot(ax12,data.freq,imag(data.s12),'DisplayName',strcat(label,' Im'))
    title(ax12,'S12')
    hold(ax12,'off')
    
    plot(ax22,data.freq,real(data.s22),'DisplayName',strcat(label,' Re'))
    plot(ax22,data.freq,imag(data.s22),'DisplayName',strcat(label,' Im'))
    title(ax22,'S22')
    hold(ax22,'off')
    
    for i=1:4
        ax = axes{i};
        grid(ax,'on');
        xlabel(ax,'Frequency (Hz)');
        legend(ax,'Orientation','horizontal','Location','southoutside');
    end
    
end