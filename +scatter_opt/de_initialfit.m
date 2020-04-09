function [p_vec,z_vec,C_vec,rss_vec] = de_initialfit(nrw,varargin)
% Perform initial fit of dielectric model to extracted eps. Poles 
% must be real, but zeros and C are allowed to be complex. Poles are
% fitted to eps with the constraint that the real and imaginary parts of
% the fitted eps cannot be negative.
% Parameters
% ----------
% nrw : table with freq and extracted mu and eps
% NPoles : Number of poles to be distributed uniformly acress frequency 
%   range. Defaults to 2
% PoleFreq: pole frequencies. If specified, NPoles is ignored, and one pole
%   is initialized at each frequency specified.
	import scatter_opt.*
    parser = inputParser;
    addOptional(parser,'NPoles',2)
    addOptional(parser,'PoleFreq',[])
    parse(parser,varargin{:})
    
    % use a subset of points for linear fit
    fit_data = nrw(1:10:801,:);
    fit_data.omega = fit_data.freq*2*pi;
    
    
    if isempty(parser.Results.PoleFreq)
        % distribute starting poles uniformly
        n_poles = parser.Results.NPoles;
        f_range = max(nrw.freq) - min(nrw.freq);
        f_inc = f_range/(n_poles+1);
        w_poles = zeros(n_poles,1);
        for n=1:n_poles
            w_poles(n) = 2*pi*(min(nrw.freq) + f_inc*n);
        end
    else
        f_poles = parser.Results.PoleFreq;
        n_poles = length(f_poles);
        w_poles = 2*pi*f_poles;
    end
    
    p_vec = zeros(1,n_poles);
    z_vec = zeros(1,n_poles);
    C_vec = zeros(1,n_poles);
    rss_vec = zeros(1,n_poles);
    y = fit_data.eps;
    
%     function rss = objfun(x,y,freq)
%         y_hat = de_eval(x,freq);
%         rss = sum(abs(y-y_hat));
%     end
    
    function cost = objfun(x,p,y,freq)
        x_full = [p x];
        y_hat = scatter_opt.de_eval(x_full,freq);
        rss = sum(abs(y-y_hat)).^2;
        % force real and imag parts to be non-negative
        % yhr = real(y_hat); yhi = imag(y_hat);
        % penalty = sum(yhr(yhr<0).^2) + sum(yhi(yhi<0).^2);
        cost = rss; %+ 1e6*penalty;
    end

    % estimate C and z for each pole
    for n=1:n_poles
        w_pole = w_poles(n);
        p = -w_pole;
        
        % initialize z at similar magnitude; set real part 100x less than
        % imag part for stability (?)
        z0 = p/100 - 1i*p;
        % estimate C
        x_pre = de_flatten_terms(p,z0,1);
        y_pred0 = de_eval(x_pre,fit_data.freq);
        C0 = median(y./y_pred0);
%         disp(mean(y./y_pred0))
%         disp(mean(abs(y./y_pred0)))
        x0 = de_flatten_terms(p,z0,C0);
        % exclude p from optimization
        x0 = x0(2:end);
        
        fun = @(x)objfun(x,p,y,fit_data.freq);
        options = optimset('MaxFunEvals',3000);
        [x_opt,fval] = fminsearch(fun,x0,options);
%         disp(fval)
        x_opt = [p x_opt];
        [p,z,C] = de_expand_x(x_opt);
        
        % linear fit estimate not good - traps in bad local min
        %-------------------------------------------------------
%         % fit real
%         yr = real(y);
%         % fixed p fixes coefficient for omega*epsilon''
%         yr = yr + (1/p)*fit_data.omega.*imag(fit_data.eps);
%         % fit remaining coefficient and intercept
%         xr = fit_data.omega;
%         fitr = fitlm(xr,yr);
%         % unpack coef and intercept from fit for convenience
%         br = fitr.Coefficients.Estimate('(Intercept)');
%         ci = fitr.Coefficients.Estimate('x1')*p;
%         
% 
%         % fit imag
%         yi = imag(y);
%         % fixed p fixes coefficient for omega*epsilon'
%         yi = yi - (1/p)*fit_data.omega.*real(fit_data.eps);
%         % fit remaining coefficient and intercept
%         xi = fit_data.omega;
%         fiti = fitlm(xi,yi);
%         % unpack coef and intercept from fit for convenience
%         cr = -fiti.Coefficients.Estimate('x1')*p;
%         bi = fiti.Coefficients.Estimate('(Intercept)');
%         
%         % extract z and C from coefficients and intercepts
%         zi = ((cr/ci)*bi*p - br*p)/(cr.^2/ci + ci);
%         zr = (bi*p-cr*zi)/ci;
% %         if abs(zr) >= abs(zi)/10
% %             zr = zr/10;
% %             cr = cr*10;
% %         end
%         C = cr + 1i*ci;
%         z = zr + 1i*zi;
        
        % grid search unhelpful - best result is always fitlm estimate
        %--------------------------------------------------------------
%         % 5-D grid search near estimated C and z: 3^4=81, 4^4=256, 5^4=625
%         p_grid = [p/2 p p*2];
%         zr_grid = [zr/10 zr/5 zr zr*5 zr*10];
%         zi_grid = [zi/10 zi/5 zi zi*5 zi*10];
%         cr_grid = [cr/10 cr/5 cr cr*5 cr*10];
%         ci_grid = [ci/10 ci/5 ci ci*5 ci*10];
%         
%         [p_grid,zr_grid,zi_grid,cr_grid,ci_grid] = ndgrid(p_grid,zr_grid,zi_grid,cr_grid,ci_grid);
%         sz = size(p_grid);
%         num = sz(1)*sz(2)*sz(3)*sz(4)*sz(5);
%         % flatten arrays
%         p_col = reshape(p_grid,[num,1]);
%         zr_col = reshape(zr_grid,[num,1]);
%         zi_col = reshape(zi_grid,[num,1]);
%         cr_col = reshape(cr_grid,[num,1]);
%         ci_col = reshape(ci_grid,[num,1]);
%         z_col = zr_col + 1i*zi_col;
%         c_col = cr_col + 1i*ci_col;
%         
%         x_cell = arrayfun(@de_flatten_terms,p_col,z_col,c_col,'Uniform',0);
%         fun = @(x)objfun(x,y,fit_data.freq);
%         f_col = cellfun(fun,x_cell);

        % lsq optimization from fitlm estimate not helpful - already at
        % local min
        %--------------------------------------------------------------
%         x0 = de_flatten_terms(p,z,C);
%         fun = @de_eval;
%         x_opt = lsqcurvefit(fun,x0,fit_data.freq,y);
        
        % store results for this pole
        p_vec(n) = p; z_vec(n) = z; C_vec(n) = C;
        
        % calculate rss for this pole
%         x = de_flatten_terms(p,z,C);
        y_pred = de_eval(x_opt,fit_data.freq);
%         disp(y_pred)
        rss_vec(n) = sum(abs(y-y_pred)).^2;
        y = y./y_pred;
    end

end
    