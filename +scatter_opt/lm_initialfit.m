function result = lm_initialfit(data,varargin)
% Perform initial fit of Laurent model with 1st-order poles via the
% rationalfit function. If NumPoles is specified, fit the specified number 
% of poles. Otherwise, perform fits with 1 to MaxNumpoles and return
% metrics for each fit to aid in determining the appropriate number of 
% poles.
% 
% Parameters
% ----------
% data : table with freq, mu, and eps
%
% Named parameters
% ----------------
% MaxNumPoles : maximum number of poles to fit
% CrossValidate : if true, perform cross-validation for each number of 
%   poles and return CV metrics
% TrainSize : fraction of data to use for training in CV. Defaults to 0.1
% EvalSplits : number of CV splits to evaluate. Default 5
% NPoles : number of poles to fit. If specified, MaxNumPoles will be 
%   ignored and fit(s) will be performed only with NumPoles
    import scatter_opt.*
    parser = inputParser;
    addParameter(parser,'MaxNumPoles',10)
    addParameter(parser,'CrossValidate',true)
    addParameter(parser,'TrainSize',0.1)
    addParameter(parser,'EvalSplits',5)
    addParameter(parser,'NPoles',-1)
    parse(parser,varargin{:})
    
    % suppress warnings from rationalfit
    warning('off','rf:rationalfit:ErrorToleranceNotMet')
    warning('off','rf:rationalfit:CheckYourData')
    
    max_num_poles = parser.Results.MaxNumPoles;
    cross_validate = parser.Results.CrossValidate;
    train_size = parser.Results.TrainSize;
    eval_splits = parser.Results.EvalSplits;
    num_poles = parser.Results.NPoles;
    
    if cross_validate 
        % set up random order for cross-validation
        num = height(data);
        rand_index = randperm(num);
        train_num = round(num*train_size);
    end
    
    if num_poles~=-1
        min_num_poles = num_poles;
        max_num_poles = num_poles;
    else
        min_num_poles = 1;
    end
    
    polrng = 1 + max_num_poles - min_num_poles;
    
    mu_cv_mean = zeros(polrng,1);
    mu_cv_min = zeros(polrng,1);
    mu_cv_std = zeros(polrng,1);
    eps_cv_mean = zeros(polrng,1);
    eps_cv_min = zeros(polrng,1);
    eps_cv_std = zeros(polrng,1);
    mu_errs = zeros(polrng,1);
    mu_scores = zeros(polrng,1);
    eps_errs = zeros(polrng,1);
    eps_scores = zeros(polrng,1);
    mu_fits = rfmodel.rational.empty(polrng,0);
    eps_fits = rfmodel.rational.empty(polrng,0);
    
    for np = min_num_poles:max_num_poles
        % Need index for cases where min_num_poles > 1
        idx = 1 + np - min_num_poles;
        if cross_validate
            mu_cv_score = double.empty;
            eps_cv_score = double.empty;
            for i = 1:eval_splits
                % get training and test data
                train_idx = rand_index((i-1)*train_num + 1:i*train_num);
                if i==1
                    % train data is first split
                    test_idx = rand_index(i*train_num + 1:num);
                elseif i==round(num*train_size)
                    % train data is last split
                    test_idx = rand_index(1:(i-1)*train_num - 1);
                else
                    % train data is in the middle
                    test_idx = [rand_index(1:(i-1)*train_num),...
                        rand_index(i*train_num + 1:num)];
                end
                
                train_data = data(train_idx,:);
                train_data = sortrows(train_data,'freq');
                
                test_data = data(test_idx,:);
                test_data = sortrows(test_data,'freq');
                
                % fit training data
                [mu_cv_fit,mu_train_score] = rationalfit(...
                    train_data.freq,train_data.mu,'NPoles',np);
                [eps_cv_fit,eps_train_score] = rationalfit(...
                    train_data.freq,train_data.eps,'NPoles',np);
                % score on test data
                mu_cv_score(i)= score_ratfit(mu_cv_fit,test_data.mu,...
                    test_data.freq);
                eps_cv_score(i) = score_ratfit(eps_cv_fit,test_data.eps,...
                    test_data.freq);
            end
            % aggregate cv scores across training splits
            mu_cv_mean(idx) = mean(mu_cv_score);
            mu_cv_std(idx) = std(mu_cv_score);
            mu_cv_min(idx) = min(mu_cv_score);
            eps_cv_mean(idx) = mean(eps_cv_score);
            eps_cv_std(idx) = std(eps_cv_score);
            eps_cv_min(idx) = min(eps_cv_score);
        end     
        % fit full dataset
        [mu_fit,mu_err] = rationalfit(data.freq,data.mu,'NPoles',np);
        [eps_fit,eps_err] = rationalfit(data.freq,data.eps,'NPoles',np);
        mu_errs(idx) = mu_err;
        mu_scores(idx) = score_ratfit(mu_fit,data.mu,data.freq);
        eps_errs(idx) = eps_err;
        eps_scores(idx) = score_ratfit(eps_fit,data.eps,data.freq);
        mu_fits(idx) = mu_fit;
        eps_fits(idx) = eps_fit;
    end
    
    % turn warnings back on
    warning('on','rf:rationalfit:ErrorToleranceNotMet')
    warning('on','rf:rationalfit:CheckYourData')
    
    if cross_validate
        result = table((min_num_poles:max_num_poles).',mu_fits.',mu_errs,...
            mu_scores,mu_cv_mean,mu_cv_min,mu_cv_std,eps_fits.',...
            eps_errs, eps_scores, eps_cv_mean,eps_cv_min,eps_cv_std,...
            'VariableNames',{'num_poles' 'mu_fit' 'mu_err'...
            'mu_score' 'mu_cv_mean' 'mu_cv_min' 'mu_cv_std'  'eps_fit'...
            'eps_err' 'eps_score' 'eps_cv_mean' 'eps_cv_min' 'eps_cv_std'});
    else        
        result = table((min_num_poles:max_num_poles).',mu_fits.',mu_errs,...
            eps_fits.',eps_errs,...
            'VariableNames',{'num_poles' 'mu_fit'...
            'mu_err' 'eps_fit' 'eps_err'});
    end
end