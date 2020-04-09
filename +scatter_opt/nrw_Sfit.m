function [s11,s21] = nrw_Sfit(Sdata,MinR2)
% fit rational model to s11 and s21
    import scatter_opt.*
	done_11 = false;
    done_21 = false;
    % suppress warnings from rationalfit
    warning('off','rf:rationalfit:ErrorToleranceNotMet')
    warning('off','rf:rationalfit:CheckYourData')
    for np=1:10
        if ~done_11
            fit11 = rationalfit(Sdata.freq,Sdata.s11,'NPoles',np);
            r2_11 = score_ratfit(fit11,Sdata.s11,Sdata.freq);
            if r2_11 > MinR2
                disp(['S11 fit: ',num2str(np),' poles, r2=',num2str(round(r2_11,5))])
                done_11 = true;
            end
        end
        if ~done_21
            fit21 = rationalfit(Sdata.freq,Sdata.s21,'NPoles',np);
            r2_21 = score_ratfit(fit21,Sdata.s21,Sdata.freq);
            if r2_21 > MinR2
                disp(['S11 fit: ',num2str(np),' poles, r2=',num2str(round(r2_21,5))])
                done_21 = true;
            end
        end
        if done_11 && done_21
            break
        end
    end 
    % turn warnings back on
    warning('on','rf:rationalfit:ErrorToleranceNotMet')
    warning('on','rf:rationalfit:CheckYourData')
    
    s11 = freqresp(fit11,Sdata.freq);
    s21 = freqresp(fit21,Sdata.freq);
    figure; ax1 = subplot(2,2,1); ax2 = subplot(2,2,2);
    ax3 = subplot(2,2,3); ax4 = subplot(2,2,4);
    axes11 = {ax1 ax2}; axes21 = {ax3 ax4};
    suptitle('S11 and S21 fits')
    title(ax1,'Re(S11)'); title(ax2,'Im(S11)')
    title(ax3,'Re(S21)'); title(ax4,'Im(S21)')
    plot_ratfit(fit11,Sdata.s11,Sdata.freq,'axes',axes11);
    plot_ratfit(fit21,Sdata.s21,Sdata.freq,'axes',axes21);
end