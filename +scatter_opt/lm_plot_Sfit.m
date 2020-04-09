function lm_plot_Sfit(lmfit,data)
    import scatter_opt.*
	pred = lm_rev_eval(lmfit,data.freq);
    
    ax11 = subplot(2,2,1);
    hold on
    plot(ax11,data.freq,real(data.s11),'.b','DisplayName','Real Measured')
    plot(ax11,pred.freq,real(pred.s11),'-r','DisplayName','Real Fit')
    plot(ax11,data.freq,imag(data.s11),'.g','DisplayName','Imag Measured')
    plot(ax11,pred.freq,imag(pred.s11),'--r','DisplayName','Imag Fit')
    title(ax11,'S11')
    hold off
    
    ax21 = subplot(2,2,2);
    hold on
    plot(ax21,data.freq,real(data.s21),'.b','DisplayName','Real Measured')
    plot(ax21,pred.freq,real(pred.s21),'-r','DisplayName','Real Fit')
    plot(ax21,data.freq,imag(data.s21),'.g','DisplayName','Imag Measured')
    plot(ax21,pred.freq,imag(pred.s21),'--r','DisplayName','Imag Fit')
    title(ax21,'S21')
    hold off
    
    ax12 = subplot(2,2,3);
    hold on
    plot(ax12,data.freq,real(data.s12),'.b','DisplayName','Real Measured')
    plot(ax12,pred.freq,real(pred.s12),'-r','DisplayName','Real Fit')
    plot(ax12,data.freq,imag(data.s12),'.g','DisplayName','Imag Measured')
    plot(ax12,pred.freq,imag(pred.s12),'--r','DisplayName','Imag Fit')
    title(ax12,'S12')
    hold off
    
    ax22 = subplot(2,2,4);
    hold on
    plot(ax22,data.freq,real(data.s22),'.b','DisplayName','Real Measured')
    plot(ax22,pred.freq,real(pred.s22),'-r','DisplayName','Real Fit')
    plot(ax22,data.freq,imag(data.s22),'.g','DisplayName','Imag Measured')
    plot(ax22,pred.freq,imag(pred.s22),'--r','DisplayName','Imag Fit')
    title(ax22,'S22')
    hold off
    
    axes = {ax11 ax21 ax12 ax22};
    for i=1:4
        ax = axes{i};
        grid(ax,'on')
        xlabel(ax,'Frequency (Hz)')
    end
    legend(ax22,'Orientation','horizontal','Location','southoutside')
end