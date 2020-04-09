function basebranch = nrw_branchSelect(Sdata,L,lambda_c,method,branchlims,...
    showplot)
% Select correct branch of permittivity and permeability
% Parameters
% ----------
% Sdata: table with freq, s11, and s21
% L: sample thickness
% lambda_c: cutoff wavelength
% method: method for branch selection. Options are GroupDelay, KK, and
%   rationalfit
% branchlims: limits of branches to consider for GroupDelay and rationalfit
%   methods
% showplot: if true, show a diagnostic plot for the branch selection
%   (differs by method)
	import scatter_opt.*
    rc = nrw_ref_coef(Sdata.s11,Sdata.s21);
    tc = nrw_trans_coef(Sdata.s11,Sdata.s21);

    if strcmp(method,'GroupDelay')
        % Group delay method. Cannot distinguish branches at high
        % frequencies
        
        % measured group delay
        tc_phase = unwrap(angle(tc));
        tmeas = (1/(2*pi))*gradient(tc_phase,Sdata.freq);
            
        branches = branchlims(1):branchlims(2);
        branch_error = zeros(length(branches),1);
        for i = 1:length(branches)
            branch = branches(i);
            lam = nrw_lambda(tc,L,branch);
            c_vac = 299792458.000176;
            lambda_0 = c_vac./Sdata.freq;
            mu = (1+rc)./(lam.*(1-rc).*sqrt((1./lambda_0.^2)-(1./lambda_c.^2)));
            eps = (lambda_0.^2./mu).*((1./lambda_c.^2)+1./lam.^2);
            % calculated group delay for branch
            calfun = sqrt(eps.*mu./(lambda_0.^2) - 1/lambda_c.^2);
            tcalc = -L*gradient(calfun,Sdata.freq);
            % compare measured and calculated group delays for lowest
            % frequencies
            num_freq = round(height(Sdata)/4);
            tdiff = tcalc - tmeas;
            % convert errors to fractions of mean tmeas (don't use
            % point-by-point relative error since tmeas oscillates about 0)
            terr = abs(tdiff)./mean(abs(tmeas));
            branch_error(i) = sum(terr(1:num_freq))/num_freq;
        end
        
        % choose the branch with the smallest error
        [min_err,min_idx] = min(branch_error);
        % check the delta between best and next best branches
        err_sort = sort(branch_error);
        next_err = err_sort(2);
        if next_err - min_err <= 0.01
            warning('Lowest error differs from next lowest error by less than 1% of tmeas. Branch choice may be ambiguous')
        end
        basebranch = branches(min_idx);
        printtable = array2table([branches.' branch_error],...
            'VariableNames',{'Branch' 'Mean_error'});
        disp(['Selected branch ', num2str(basebranch),...
            ' with mean group delay error ', num2str(min_err)])
        disp('All branch group delay errors:')
        disp(printtable)
        
        if showplot
            % plot tmeas and tcalc for chosen branch
            lam = nrw_lambda(tc,L,basebranch);
            mu = (1+rc)./(lam.*(1-rc).*sqrt((1./lambda_0.^2)-(1./lambda_c.^2)));
            eps = (lambda_0.^2./mu).*((1./lambda_c.^2)+1./lam.^2);
            % calculated group delay for branch
            calfun = sqrt(eps.*mu./(lambda_0.^2) - 1/lambda_c.^2);
            tcalc = -L*gradient(calfun,Sdata.freq);
            figure; axes = plot_rivf(Sdata.freq,tmeas,'label',...
                'Meas.');
            axes = plot_rivf(Sdata.freq,tcalc,axes,'label',...
                'Calc.');
            suptitle('Measured vs. Calculated Group Delay')
            legend(axes{1});
        end
            
    elseif strcmp(method,'KK')
        % Kramers-Kronig method using refractive index. Not recommended for
        % finite frequency ranges of small width
        % calculate refractive index for branch 0
        n_eff = nrw_refracIndex(Sdata,L);
        % get imag part of n_eff - no branch ambiguity
        kap = imag(n_eff);
        % construct real part from imag part based on Kramers-Kronig
        n_re_KK = nrw_refracKK(Sdata.freq,kap);
        % free space wavenumber
        c_vac = 299792458.000176;
        k0 = 2*Sdata.freq./c_vac;
        % correct branch given by branch with n_re closest to n_re_KK
        kk_branch = round((n_re_KK - real(n_eff)).*k0*L/(2*pi));
        % Use the KL to find the initial branch based on most common branch
        % at low frequency
        num_freq = round(height(Sdata)/8);
        % avoid lowest frequencies - truncation errors
        basebranch = mode(kk_branch(10:num_freq));
        disp(['Selected branch ',num2str(basebranch)])
        
        if showplot
            % plot n_re_KK vs. n_re for basebranch +/- 1
            figure; plot(Sdata.freq,n_re_KK,'DisplayName','KK');
            hold on
            for b = basebranch-1:basebranch+1
                plot(Sdata.freq,real(n_eff) + 2*b*pi./(k0*L),...
                    'DisplayName',['Branch ' num2str(b)]);
            end
            hold off
            legend
            title('Refractive index: KK vs. branches')
        end
            
    elseif strcmp(method,'rationalfit')
        alpha = -log(abs(tc))./L;
        beta0 = unwrap((-angle(tc)))/L;
        num_poles = 1:5;
        branches = branchlims(1):branchlims(2);
        gam_fits = array2table(num_poles.','VariableNames',{'NPoles'});
        
        % suppress warnings from rationalfit
        warning('off','rf:rationalfit:ErrorToleranceNotMet')
        warning('off','rf:rationalfit:CheckYourData')
    
        branch_med_r2 = zeros(length(branches),1);
        for i = 1:length(branches)
            branch = branches(i);
            if branch < 0
                bname = strcat('bn',num2str(abs(branch)));
            else
                bname = strcat('b',num2str(branch));
            end
            gam_fits.(strcat(bname,'_err')) = zeros(height(gam_fits),1);
            gam_fits.(strcat(bname,'_ReR2')) = zeros(height(gam_fits),1);
            gam_fits.(strcat(bname,'_ImR2')) = zeros(height(gam_fits),1);
            gam_fits.(strcat(bname,'_MeanR2')) = zeros(height(gam_fits),1);
            gam_fits.(strcat(bname,'_fit')) = rfmodel.rational.empty(length(num_poles),0);
            
            for np = num_poles
                gamma = alpha + 1i*(beta0 + 2*pi*branch/L);
                [fit,err] = rationalfit(Sdata.freq,gamma,'NPoles',np);
                gam_fits.(strcat(bname,'_err'))(np) = err;
                gam_fits.(strcat(bname,'_fit'))(np) = fit;
                re_r2 = score_ratfit(fit,gamma,Sdata.freq,'Part','real');
                im_r2 = score_ratfit(fit,gamma,Sdata.freq,'Part','imag');
                gam_fits.(strcat(bname,'_ReR2'))(np) = re_r2;
                gam_fits.(strcat(bname,'_ImR2'))(np) = im_r2;
                gam_fits.(strcat(bname,'_MeanR2'))(np) = (re_r2 + im_r2)/2;
            end
            branch_med_r2(i) = median(gam_fits.(strcat(bname,'_MeanR2')));
        end
        disp(gam_fits)
        % choose the branch that maximize the median r2 score
        [r2max,max_idx] = max(branch_med_r2);
        r2_sort = sort(branch_med_r2);
        r2_next = r2_sort(end-1);
        if r2max - r2_next <= 0.05
            warning('Best R2 score is less than 0.05 greater than next best. Branch choice may be ambigous')
        end
        basebranch = branches(max_idx);
        printtable = array2table([branches.' branch_med_r2],...
            'VariableNames',{'Branch' 'Median_R2'});
        disp(['Selected branch ', num2str(basebranch),...
            ' with median R2 score ', num2str(r2max)])
        disp('All branch R2 scores:')
        disp(printtable)
        
        % turn warnings back on
        warning('on','rf:rationalfit:ErrorToleranceNotMet')
        warning('on','rf:rationalfit:CheckYourData')
        
        if showplot
            % plot median rationalfit of chosen branch +/- 1
            if basebranch < 0
                bbname = strcat('bn',num2str(abs(basebranch)));
            else
                bbname = strcat('b',num2str(basebranch));
            end
            row_idx = find(gam_fits.(strcat(bbname,'_MeanR2'))==r2max,1);
            figure;
            axr = subplot(1,2,1); hold(axr,'on')
            axi = subplot(1,2,2); hold(axi,'on')
            
            % real part is same for all
            plot(axr,Sdata.freq,alpha,'k','DisplayName','(all)');
            
            colors = ['b' 'r' 'g'];
            i = 1;
            for b = max([basebranch-1 min(branches)]):min([basebranch+1 max(branches)])
                if b < 0
                    bname = strcat('bn',num2str(abs(b)));
                else
                    bname = strcat('b',num2str(b));
                end
                fit = gam_fits.(strcat(bname,'_fit'))(row_idx);
                gamma = alpha + 1i*(beta0 + 2*pi*b/L);
                [pred,outfreq] = freqresp(fit,Sdata.freq);
                c = colors(i);
                
                plot(axi,Sdata.freq,imag(gamma),c,'DisplayName',...
                    ['Br. ' num2str(b)]);
                plot(axr,Sdata.freq,real(pred),[c '--'],'DisplayName',...
                    ['Br. ' num2str(b)]);
                plot(axi,Sdata.freq,imag(pred),[c '--'],...
                    'HandleVisibility','off')

                i = i+1;
            end
            legend(axr)
            legend(axi)
            xlabel(axr,'Frequency(HZ)')
            xlabel(axi,'Frequency(HZ)')
            ylabel(axr,'\alpha')
            ylabel(axi,'\beta')
            suptitle('Rational fits of \gamma')
        end

    else
        error('Invalid branch selection method. Options are GroupDelay, KK, and rationalfit')
    end