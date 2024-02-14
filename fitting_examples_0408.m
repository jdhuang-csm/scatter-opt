% 40NZF
nzf4 = Sproc(1.64e-3,1.6e-3,0.04572);
nzf4.load('data/X_40NZF_notape_Ag_31Jan.txt','HeaderLines',8)
% plot measured S params
nzf4.plot_Sij
% extract using rationalfit branch selection, and plot the branch selection
% results
nzf4.extract('nrw','PlotBranchSelect',true)
% extract with optimization of L1 and L2
nzf4.extract('nrw','PlotBranchSelect',true,'OptimizeL',true)
% extract after fitting rational model to S parameters
nzf4.extract('nrw','PlotBranchSelect',true,'FitS',true)
% extract with manual branch selection
nzf4.extract('nrw','BaseBranch',0)

% for fit of 40NZF - use FitS
nzf4.extract('nrw','PlotBranchSelect',true,'FitS',true)

% plot extracted PP
nzf4.plot_PP
% plot S params calculated from extracted PP - should match exactly if
% using NRW without OptimizeL or FitS
nzf4.plot_Sreverse

% fit LM model to extracted PP for 1-10 poles
nzf4.lm_initialfit

% Get x0 using 7 poles each for mu and eps
nzf4.lm_estimate_x0('NPoles1',7)
% Get x0 using 7 poles each for mu and eps, without estimating 2nd-order
% poles (generally this results in the 2nd-order poles having no effect)
nzf4.lm_estimate_x0('NPoles1',7,'EstOrd2',false)

% fit LM model
nzf4.lm_lsqfit
% fit LM model with shifted 1st-order poles
nzf4.lm_lsqfit('ShiftFOPoles',true)
% fit LM model with more FEvals or iterations allowed
nzf4.lm_lsqfit('MaxFEval',5000,'MaxIter',350)

% plot fit results
nzf4.lm_plot_Sfit
nzf4.lm_plot_PPfit