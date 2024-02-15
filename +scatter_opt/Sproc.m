classdef Sproc < handle
       properties
           L {mustBeNumeric}
           L_air {mustBeNumeric}
           lambda_c {mustBeNumeric}
           Sdata % converted S-parameter data
           PPdata % extracted mu and epsilon
           extract_method % method used for extraction; read-only
           mu_np1 % number of 1st-order poles for mu
           eps_np1 % number of 1st-order poles for eps
           np2 % number of 2nd-order poles
           wg_x0 % waveguide initial parameters
           mu_x0 % mu initial parameters
           eps_x0 % eps initial parameters
           lm_initfits %
           lmfit % lm fit object
           lm_internal % lm fit internal troubleshooting object
           defit % dielectric fit object
           de_internal % dielectric fit internal troubleshooting object
           
       end
       methods
           function obj = Sproc(L,L_air,lambda_c)
               obj.L = L;
               obj.L_air = L_air;
               obj.lambda_c = lambda_c;
           end
       % Data processing
           function load(obj,file,varargin)
			   import scatter_opt.*
               parser = inputParser;
               addOptional(parser,'HeaderLines',0)
               addOptional(parser,'Delimiter',' ')
               addParameter(parser,'PhaseCorrect',false)
               parse(parser,varargin{:})
               obj.Sdata = loadconvert_raw(file,...
                   'HeaderLines',parser.Results.HeaderLines,...
                   'Delimiter',parser.Results.Delimiter);
               if parser.Results.PhaseCorrect
                   % apply phase correction to S21 and S12 to account for 
                   % transmission line length displaced by sample
                   gam0 = rev_gamma_0(obj.Sdata.freq*2*pi,obj.lambda_c);
                   obj.Sdata.s21 = obj.Sdata.s21.*exp(-gam0*obj.L);
                   obj.Sdata.s12 = obj.Sdata.s12.*exp(-gam0*obj.L);
               end
           end
           function extract(obj,meth,varargin)
           % Extract mu and eps from S parameters using NRW or NNI method
           % Parameters
           % ----------
           % meth: 'nrw' or 'nni'
           %
           % Named (optional) parameters
           % ---------------------------
           % BaseBranch: branch of complex root to start from. Defaults to 0, will 
           %   raise warning if principal branch is verfiably incorrect
           % IncrementBranch: if true, check for discontinuities indicative of branch
           %   changes and increment branch from basebranch accordingly
				import scatter_opt.*
                if meth=='nrw'
                    obj.PPdata = nrw_extract(obj.Sdata,obj.L,obj.L_air,...
                        obj.lambda_c,varargin{:});
                elseif meth=='nni'
                    obj.PPdata = nni_extract(obj.Sdata,obj.L,obj.L_air,...
                        obj.lambda_c,varargin{:});
                else
                    error('method must be nrw or nni')
                end
                obj.extract_method = meth;
           end
       % Data plotting
           function plot_PP(obj)
           % Plot extracted mu and epsilon vs. frequency
			   import scatter_opt.*
               figure;plot_PP(obj.PPdata);
           end
           function plot_Sij(obj)
           % Plot all measured S parameters vs. frequency
			   import scatter_opt.*
               figure;plot_Sij(obj.Sdata);
           end
           function plot_Sreverse(obj)
           % plot the S parameters calculated from the reverse transform 
           % of the extracted mu and eps alongside the measured S parameters
			   import scatter_opt.*
               rev_out = rev_transform(obj.PPdata,obj.lambda_c,obj.L,0,0);
               figure; ax1 = subplot(2,2,1); ax2 = subplot(2,2,2);
               ax3 = subplot(2,2,3); ax4 = subplot(2,2,4);
               axes11 = {ax1 ax2}; axes21 = {ax3 ax4};
               plot_rivf(obj.Sdata.freq,obj.Sdata.s11,axes11,'label','Measured');
               plot_rivf(rev_out.freq,rev_out.s11,axes11,'label','Reverse');
               legend(ax1)
               title(ax1,'S11 Real')
               title(ax2,'S11 Imag')
               plot_rivf(obj.Sdata.freq,obj.Sdata.s21,axes21,'label','Measured');
               plot_rivf(rev_out.freq,rev_out.s21,axes21,'label','Reverse');
               legend(ax3)
               title(ax3,'S21 Real')
               title(ax4,'S21 Imag')
           end
       % Laurent model fitting
           function lm_initialfit(obj,varargin)
           % Perform initial fit of Laurent model with 1st-order poles via the
           % rationalfit function. If NumPoles is specified, fit the specified number 
           % of poles. Otherwise, perform fits with 1 to MaxNumpoles and return
           % metrics for each fit to aid in determining the appropriate number of 
           % poles.
           %
           % Named (optional) parameters
           % ----------------
           % MaxNumPoles : maximum number of poles to fit
           % CrossValidate : if true, perform cross-validation for each number of 
           %   poles and return CV metrics
           % TrainSize : fraction of data to use for training in CV. Defaults to 0.1
           % EvalSplits : number of CV splits to evaluate. Default 5
           % NPoles : number of poles to fit. If specified, MaxNumPoles will be 
           %   ignored and fit(s) will be performed only with NumPoles
               import scatter_opt.*
			   initfits = lm_initialfit(obj.PPdata,varargin{:})
               obj.lm_initfits = initfits;
               
           end
           function lm_plot_initfit(obj,field,NPoles)
           % Plot initial first-order fit of mu or epsilon
           % Parameters
           % ----------
           % field: which field to plot ('mu' or 'eps')
           % NPoles: number of poles
               import scatter_opt.*
			   if strcmp(field,'mu')
                   figure;plot_ratfit(obj.lm_initfits.mu_fit(NPoles),...
                       obj.PPdata.mu,obj.PPdata.freq,'FlipIm',true)
               elseif strcmp(field,'eps')
                   figure;plot_ratfit(obj.lm_initfits.eps_fit(NPoles),...
                       obj.PPdata.eps,obj.PPdata.freq,'FlipIm',true)
               else
                   error('field must be ''mu'' or ''eps''')
               end
           end
           function lm_estimate_x0(obj,varargin)
           % Estimate initial parameter vectors for mu and eps including 
           % first- and second-order poles.
           % Named (optional) parameters
           % ---------------------------
           % NPoles1: number of 1st-order poles to use for both mu and eps
           % NPoles2: number of 2nd-order poles to use for both mu and eps
           % EstOrd2: if true, estimate initial parameters for 2nd-order poles. If
           %   false, initialize 2nd-order poles parameters at zero
           % MuNPoles1: number of 1st-order poles to use for mu. If specified,
           %   overrides NPoles1
           % EpsNPoles1: number of 1st-order poles to use for eps. If specified,
           %   overrides NPoles1
			   import scatter_opt.*
               [obj.mu_x0,obj.eps_x0,obj.mu_np1,obj.eps_np1,obj.np2] = ...
                   lm_estimate_x0(obj.PPdata,varargin{:});
           end 
           
           function lm_plot_initSfit(obj)
           % Plot initial fit of S parameters
                import scatter_opt.*
				tmpfit = struct;
                tmpfit.x_mu = obj.mu_x0;
                tmpfit.x_eps = obj.eps_x0;
                tmpfit.mu_np1 = obj.mu_np1;
                tmpfit.eps_np1 = obj.eps_np1;
                tmpfit.mu_np2 = obj.np2;
                tmpfit.eps_np2 = obj.np2;
                tmpfit.lambda_c = obj.lambda_c;
                tmpfit.L = obj.L; tmpfit.L1 = 0; tmpfit.L2 = obj.L_air-obj.L;
                figure; lm_plot_Sfit(tmpfit,obj.Sdata)
           end
           
           function lm_plot_initPPfit(obj)
           % Plot initial fit of mu and epsilon
                import scatter_opt.*
				tmpfit = struct;
                tmpfit.x_mu = obj.mu_x0;
                tmpfit.x_eps = obj.eps_x0;
                tmpfit.mu_np1 = obj.mu_np1;
                tmpfit.eps_np1 = obj.eps_np1;
                tmpfit.mu_np2 = obj.np2;
                tmpfit.eps_np2 = obj.np2;
                tmpfit.lambda_c = obj.lambda_c;
                tmpfit.L = obj.L; tmpfit.L1 = 0; tmpfit.L2 = 0;
                figure; lm_plot_PPfit(tmpfit,obj.PPdata)
           end
           
           function lm_lsqfit(obj,varargin)
           % Fit LM model to S parameters
           % Named (optional) parameters
           % ----------------
           % MaxFEval : maximum number of function evaluations to allow the optimizer
           %   to perform
           % MultiStart: if true, start the optimization from multiple initial
           %   points
           % MultiStartPoints : number of points from which to initialize the
           %   optimization
           % ShiftFOPoles : if true, shift the parameters for first-order poles when 
           %   generating MultiStart points. If false, leave the first-order pole
           %   parameters the same for all MultiStart points (only shift 2nd-order
           %   poles)
           % NCCP2 : number of complex conjugate pairs in second-order poles. 0 (real 
           %   2nd-order poles only) tends to work well, but in some cases complex
           %   2nd-order poles may be beneficial.
               import scatter_opt.*
			   obj.wg_x0 = [obj.L 0 obj.lambda_c];
               [obj.lmfit,obj.lm_internal] = lm_lsqfit(obj.Sdata,obj.L_air,obj.L,...
                   obj.mu_x0,obj.eps_x0,obj.wg_x0,obj.mu_np1,obj.eps_np1,...
                   obj.np2,varargin{:});
           end
           function lm_plot_Sfit(obj)
           % Plot LM fit of S parameters
               import scatter_opt.*
			   figure;lm_plot_Sfit(obj.lmfit,obj.Sdata)
           end
           function lm_plot_PPfit(obj)
           % Plot LM fit of mu and eps
               import scatter_opt.*
			   figure;lm_plot_PPfit(obj.lmfit,obj.PPdata)
           end
           function out = lm_output(obj,varargin)
           % Output data and model results. Save to file if specified
           % Named (optional) parameters
           % ---------------------------
           % SaveFile : file to save output to
               import scatter_opt.*
			   out = lm_output(obj.lmfit,obj.Sdata,obj.PPdata,...
                   'ExtractMethod',obj.extract_method,varargin{:});
           end
               
       % Dielectric model fitting
           function de_estimate_x0(obj,varargin)
           % Perform initial fit of dielectric model to extracted eps. Poles 
           % must be real, but zeros and C are allowed to be complex. Poles are
           % fitted to eps with the constraint that the real and imaginary parts of
           % the fitted eps cannot be negative.
           % Named (optional) parameters
           % ---------------------------
           % NPoles : Number of poles to be distributed uniformly across frequency 
           %   range. Defaults to 2
           % PoleFreq: pole frequencies. If specified, NPoles is ignored, and one pole
           %   is initialized at each frequency specified.
               import scatter_opt.*
			   obj.eps_x0 = de_estimate_x0(obj.PPdata,varargin{:});
           end
           function de_plot_initPPfit(obj)
               % Plot initial fit of mu and epsilon
                import scatter_opt.*
				tmpfit = struct;
                tmpfit.x_eps = obj.eps_x0;
                figure; de_plot_PPfit(tmpfit,obj.PPdata)
           end
           function de_plot_initSfit(obj)
               % Plot initial fit of mu and epsilon
                import scatter_opt.*
				tmpfit = struct;
                tmpfit.x_eps = obj.eps_x0;
                tmpfit.lambda_c = obj.lambda_c;
                tmpfit.L = obj.L; tmpfit.L1 = 0; tmpfit.L2 = obj.L_air-obj.L;
                figure; de_plot_Sfit(tmpfit,obj.Sdata)
           end

           function de_lsqfit(obj,varargin)
           % Fit dielectric model to S parameters
           % Named (optional) parameters
           % ----------------
           % MaxFEval : maximum number of function evaluations to allow the optimizer
           %   to perform
           % MultiStart: if true, start the optimization from multiple initial
           %   points
           % MultiStartPoints : number of points from which to initialize the
           %   optimization
               import scatter_opt.*
			   obj.wg_x0 = [obj.L 0 obj.lambda_c];
               [obj.defit,obj.de_internal] = de_lsqfit(obj.Sdata,obj.L_air,...
                   obj.L,obj.eps_x0,obj.wg_x0,varargin{:});
           end
           function de_plot_Sfit(obj)
           % Plot dielectric fit of S parameters
               import scatter_opt.*
			   figure;de_plot_Sfit(obj.defit,obj.Sdata)
           end
           function de_plot_PPfit(obj)
           % Plot dielectric fit of S parameters
               import scatter_opt.*
			   figure;de_plot_PPfit(obj.defit,obj.PPdata)
           end        
           function out = de_output(obj,varargin)
           % Output data and model results. Save to file if specified
           % Named (optional) parameters
           % ---------------------------
           % SaveFile : file to save output to
               import scatter_opt.*
			   out = de_output(obj.defit,obj.Sdata,obj.PPdata,...
                   'ExtractMethod',obj.extract_method,varargin{:});
           end        
       end
end
           