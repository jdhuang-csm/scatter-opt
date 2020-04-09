function out = nni_extract(Sdata,L,lambda_c,varargin)
% Compute mu and epsilon from s11 and s21 following method in Rhode &
% Schwarz Annex 4.
% Parameters
% ----------
% data : table with freq, s11, and s21
% L : sample thickness in m
% lambda_c : cutoff wavelength
% BaseBranch: branch of complex root. Defaults to 0, will raise warning if 
%   principal branch is verfiably incorrect
% IncrementBranch: if true, check for discontinuities indicative of branch
%   changes and increment branch from basebranch accordingly
% InterpolatePoints: if true, check for remaining discontinuities after
%   incrementing branch and interpolate those points (assume that these are
%   due to branch steps being slightly offset)
	import scatter_opt.*
    parser = inputParser;
    addOptional(parser,'BranchSelectMethod','rationalfit')
    addOptional(parser,'PlotBranchSelect',false)
    addParameter(parser,'BranchLimits',[-2 2])
    addParameter(parser,'BaseBranch',0)
    addParameter(parser,'GeoMean',false)
    addParameter(parser,'OptimizeL',false)
    addParameter(parser,'FitS',false)
    addParameter(parser,'SFitMinR2',0.99)
    addParameter(parser,'PhaseUnwrap',true)
%     addParameter(parser,'IncrementBranch',true)
%     addParameter(parser,'InterpolatePoints',true)
    parse(parser,varargin{:})
    
    data = Sdata;
    
    if parser.Results.GeoMean
        % use the geometric mean of s11 and s22 in place of s11
        sii_mean = sqrt(data.s11.*data.s22);
        sign_change = sign(real(sii_mean))./sign(real(data.s11));
        data.s11 = sii_mean.*sign_change;    
    end
    
    if parser.Results.OptimizeL
        [L1,L2] = nrw_optimizeL(data,L,L_air,lambda_c);
        disp([L1,L2])
        gam0 = rev_gamma_0(2*pi*data.freq,lambda_c);
        R1 = exp(-gam0*L1); R2 = exp(-gam0*L2);
        data.s11 = mean([data.s11./R1.^2 data.s22./R2.^2],2);
        data.s21 = data.s21./(R1.*R2);
    end
    
    if parser.Results.FitS
        [s11,s21] = nrw_Sfit(data,parser.Results.SFitMinR2);
        data.s11 = s11;
        data.s21 = s21;
    end
    
    if strcmp(parser.Results.BranchSelectMethod,'manual') || any(strcmp(varargin,'BaseBranch'))
        % Use user-specified branch if BranchSelectMethod is manual OR if
        % BaseBranch was explicitly passed
        basebranch = parser.Results.BaseBranch;
        if any(strcmp(varargin,'BaseBranch')) &&... 
                any(strcmp(varargin,'BranchSelectMethod')) &&...
                ~strcmp(parser.Results.BranchSelectMethod,'manual')
            % if user specified both BaseBranch and non-manual 
            % BranchSelectMethod, raise warning
            warning('Both BranchSelectMethod and BaseBranch were specified. Assuming manual selection and using specified BaseBranch')
        end
    else
        basebranch = nrw_branchSelect(data,L,lambda_c,...
            parser.Results.BranchSelectMethod,parser.Results.BranchLimits,...
            parser.Results.PlotBranchSelect);
    end
    
    rc = nrw_ref_coef(data.s11,data.s21);
    tc = nrw_trans_coef(data.s11,data.s21);
    lam = nrw_lambda(tc,L,basebranch,parser.Results.PhaseUnwrap);
    c_vac = 299792458.000176;
    lambda_0 = c_vac./data.freq;
    lambda_og = 1./sqrt(1./lambda_0.^2 - 1./lambda_c.^2);
    n = 1;
    mu = ones(height(data),1);
    eps_eff = ((1-rc)./(1+rc)).^(n-1).*(lambda_og./lam).^(n+1);
    eps = (1-lambda_0.^2./lambda_c.^2).*eps_eff + lambda_0.^2./lambda_c.^2;
    freq = data.freq;
    
    % check for discontinuities
    ediff = diff(eps);
    step_index = find(abs(ediff) > 20*mean(abs(ediff)));
    % if contiguous sequences exist, only take the first index in each
    % contiguous sequence
    step_expand = [0; step_index];
    step_index = step_index(diff(step_expand)>1);
    if ~isempty(step_index)
        warning(['Found discontinuity(ies) in epsilon after phase unwrapping. Check data and BaseBranch. Step indexes: ',mat2str(step_index)])
    end
    
% Manual discontinuity identification replaced by phase unwrapping
% log(1/tc) in nrw_lambda
%     % check for discontinuities
%     if parser.Results.IncrementBranch
%         ediff = diff(eps);
%         step_index = find(abs(ediff) > 20*mean(abs(ediff)));
%         % if contiguous sequences exist, only take the first index in each
%         % contiguous sequence
%         step_expand = [0; step_index];
%         step_index = step_index(diff(step_expand)>1);
%         if ~isempty(step_index)
%             warning(['Found discontinuity in epsilon. Incrementing branch. Step indexes:',mat2str(step_index)])
%             % start the branch at basebranch
%             branch = zeros(height(data),1);
%             branch = branch + basebranch;
%             for j=1:length(step_index)
%                 % for each sign change, increment the branch by 1
%                 idx = step_index(j);
%                 branch(idx+1:end) = branch(idx+1:end) + 1;
%             end
%             % recalculate lambda
%             lam = nrw_lambda(tc,L,branch);
%             % recalculate eps
%             n = 1;
%             eps_eff = ((1-rc)./(1+rc)).^(n-1).*(lambda_og./lam).^(n+1);
%             eps = (1-lambda_0.^2./lambda_c.^2).*eps_eff + ...
%                 lambda_0.^2./lambda_c.^2;
%         end
%         
%         % after incrementing branch, check for remaining discontinuities
%         % (may exist due to branches stepping at slightly different spots)
%         if parser.Results.InterpolatePoints
%             ediff = diff(eps);
%             step_index = find(abs(ediff) > 20*mean(abs(ediff)));
%             if ~isempty(step_index)
%                 warning(['Found remaining discontinuity in epsilon after incrementing branch. Removing points. Step indexes:',mat2str(step_index)])
%                 % if contiguous sequences exist, find the start and end of 
%                 % each contiguous sequence
%                 step_expand = [0; step_index];
%                 start_index = find(diff(step_expand)>1);
%                 end_index = start_index - 1;
%                 end_index = end_index(2:end);
%                 end_index = [end_index; length(step_index)];
%                 % indexes are diff indexes. Increase by 1 to align with
%                 % data
%                 start_index = step_index(start_index) + 1;
%                 end_index = step_index(end_index) + 1;
%                 for j=1:length(start_index)
%                     % for each discontinuity or sequence of
%                     % discontinuities, replace the offending point(s) with
%                     % a linear interpolation of the surrounding points
%                     start_idx = start_index(j);
%                     end_idx = end_index(j);
%                     f_pre = freq(start_idx - 1);
%                     f_post = freq(end_idx + 1);
%                     mu_pre = mu(start_idx - 1);
%                     mu_post = mu(end_idx + 1);
%                     eps_pre = eps(start_idx - 1);
%                     eps_post = eps(end_idx + 1);
%                     mu_slope = (mu_post - mu_pre)/(f_post - f_pre);
%                     eps_slope = (eps_post - eps_pre)/(f_post - f_pre);
%                     mu(start_idx:end_idx) = mu_pre + ...
%                         mu_slope*(freq(start_idx:end_idx) - f_pre);
%                     eps(start_idx:end_idx) = eps_pre + ...
%                         eps_slope*(freq(start_idx:end_idx) - f_pre);
%                 end
%             end
%         end
%     end
    
    out = array2table([freq mu eps],'VariableNames',{'freq' 'mu' 'eps'});
    
end
    