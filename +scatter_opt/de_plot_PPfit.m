function de_plot_PPfit(defit,data)
% Plot fit of permittivity and permeability
% Parameters 
% ----------
% defit : dielectric model fit struct
% data : table with freq, mu, and eps
    import scatter_opt.*
    mu_pred = ones(height(data),1);
    eps_pred = de_eval(defit.x_eps,data.freq);
    
    ax11 = subplot(2,2,1);
    hold on
    plot(ax11,data.freq,real(data.mu),'-b','DisplayName','Extracted')
    plot(ax11,data.freq,real(mu_pred),'-r','DisplayName','Model')
    title(ax11,'Re(mu)')
    hold off
    
    ax12 = subplot(2,2,2);
    hold on
    plot(ax12,data.freq,-imag(data.mu),'-b','DisplayName','Extracted')
    plot(ax12,data.freq,-imag(mu_pred),'-r','DisplayName','Model')
    title(ax12,'Im(mu)')
    hold off
    
    ax21 = subplot(2,2,3);
    hold on
    plot(ax21,data.freq,real(data.eps),'-b','DisplayName','Extracted')
    plot(ax21,data.freq,real(eps_pred),'-r','DisplayName','Model')
    title(ax21,'Re(eps)')
    hold off
    
    ax22 = subplot(2,2,4);
    hold on
    plot(ax22,data.freq,-imag(data.eps),'-b','DisplayName','Extracted')
    plot(ax22,data.freq,-imag(eps_pred),'-r','DisplayName','Model')
    title(ax22,'Im(eps)')
    hold off
    
    axes = {ax11 ax12 ax21 ax22};
    for i=1:4
        ax = axes{i};
        grid(ax,'on')
        xlabel(ax,'Frequency (Hz)')
        legend(ax,'Orientation','horizontal','Location','southoutside')
    end
end