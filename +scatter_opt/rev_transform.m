function varargout = rev_transform(nrw_out,lambda_c,L,L1,L2,varargin)
% Calculate s11, s21, s12, and s22 from NRW-calculated mu and epsilon
% Parameters
% ----------
% nrw_out : output of nrw_extract
% lambda_c : cutoff wavelength
% L : sample thickness
% L1: distance between reference plane and sample face on side 1
% L2: distance between reference plane and sample face on side 2
% conv: converted raw data table with freq, s11, s21. Required only if
%   gam_sign and gam0_sign need to be determined automatically.
% gam_sign : sign of gamma (+1 for positive, -1 for negative). If 0, 
%   determine automatically. 
% gam0_sign : sign of gamma_0 (+1 for positive, -1 for negative). If 0, 
%   determine automatically.
    import scatter_opt.*
	parser = inputParser;
%     addOptional(parser,'conv',array2table(zeros(0,0)))
%     addOptional(parser,'gam_sign',0)
%     addOptional(parser,'gam0_sign',0)
    addParameter(parser,'modes',1)
    addParameter(parser,'lambda_c2',0)
    addParameter(parser,'betas',[])
    parse(parser,varargin{:})
    
%     if parser.Results.gam_sign==0
%         if height(parser.Results.conv)==0
%             error('Table conv must be supplied to determine gam_sign and gam0_sign. See documentation')
%         else
%             conv = parser.Results.conv;
%         end
%     end
    
    omega = 2.*pi.*nrw_out.freq;
    R1 = rev_Ri(omega,lambda_c,L1);
    R2 = rev_Ri(omega,lambda_c,L2);
    gam = rev_gamma(omega,nrw_out.mu,nrw_out.eps,lambda_c);
    
    gam0 = rev_gamma_0(omega,lambda_c);
%     disp(sum(sign(real(gam))))
    
    % ensure real part of gamma (attenuation constant) is positive
    gr_sign = sign(real(gam));
    % if real part is exactly zero, treat as positive
    gr_sign(gr_sign==0) = 1;
%     disp(sum(gr_sign))
    gam = gam.*gr_sign;
    
    % phase unwrap beta
%     disp(gam(330:360));
    beta = imag(gam)*L;
    buw = unwrap(beta);
    gam = real(gam) + 1i*buw/L;
%     disp(beta(330:360))
    
    
    % match sign of gam0 to sign of gam 
    % this seems to be necessary to ensure that rev_transform of the NRW
    % matches the measured S params. However, it doesn't seem to change the
    % final fit results (the same solution can be found either way)
%     gam0 = gam0.*gr_sign;
    % causes steps for teflon!

    
   
    % get correct signs for gamma and gamma_0 - real parts (attenuation 
    % constants) must be non-negative
    
%     if real(gam(1)) < 0
%         if abs(sum(sign(real(gam))))==height(nrw_out)
%             gam = -gam;
%             % apply the same sign to gam0 (this seems to work? not sure why
%             % gam0 sign would change)
%             gam0 = -gam0;
%         else
% %             warning('Gamma changes sign')
%             gam = -gam;
%             % apply the same sign to gam0 (this seems to work? not sure why
%             % gam0 sign would change)
%             gam0 = -gam0;
%         end
%     end


%     
%     if abs(sum(sign(real(gam0))))==height(nrw_out)
%         if real(gam0(1)) < 0
%             gam0 = -gam0;
%         end
%     elseif sum(abs(sign(real(gam0))))==0
%         % real part is zero (expected)
%         gam0 = gam0;
%     else
%         error('Gamma_0 changes sign')
%     end
    
    % get correct sign of gamma by matching rev_trans_coef to
    % nrw_trans_coef
%     if parser.Results.gam_sign==0
%         % correct sign may switch - sign must be a vector
%         gam_sign = zeros(length(nrw_out.freq),1);
%         disp(sum(sign(real(gam))))
%         tc_nrw = nrw_trans_coef(conv.s11,conv.s21);
%         tcp = rev_trans_coef(gam,L);
%         tcm = rev_trans_coef(-gam,L);
%         % when tcp matches tc_nrw, sign is +1
%         gam_sign(abs(tc_nrw - tcp)./abs(tc_nrw) < 1e-3) = 1;
%         % when tcm matches tc_nrw, sign is -1
%         gam_sign(abs(tc_nrw - tcm)./abs(tc_nrw) < 1e-3) = -1;
%         
% %         if max(abs(tc_nrw - tcp)./abs(tc_nrw)) < 1e-3
% %             gam_sign = 1;
% %         elseif max(abs(tc_nrw - tcm)./abs(tc_nrw)) < 1e-3
% %             gam_sign = -1;
%         nomatch = sum(gam_sign==0);
%         if nomatch > 0
%             % If trans coefs don't match with either sign, warn and perform
%             % calculation with positive sign
%             error('Reverse transmission coefficient does not match NRW transmission coefficient')
% %             gam_sign = 1;
%         end
%     else
%         gam_sign = parser.Results.gam_sign;
%     end
%     gam = gam_sign.*gam;

    tc = rev_trans_coef(gam,L);
%     disp(tc(330:360))
    
%     % get correct sign of gamma_0 by matching rev_ref_ceof to nrw_ref_coef
%     if parser.Results.gam0_sign==0
%         % correct sign may switch - sign must be a vector
%         gam0_sign = zeros(length(nrw_out.freq),1);
%         rc_nrw = nrw_ref_coef(conv.s11,conv.s21);
%         rcp = rev_ref_coef(nrw_out.mu,gam0,gam);
%         rcm = rev_ref_coef(nrw_out.mu,-gam0,gam);
%         % when rcp matches rc_nrw, sign is +1
%         gam0_sign(abs(rc_nrw - rcp)./abs(rc_nrw) < 1e-3) = 1;
%         % when rcm matches rc_nrw, sign is -1
%         gam0_sign(abs(rc_nrw - rcm)./abs(rc_nrw) < 1e-3) = -1;
% %         if max(abs(rc_nrw - rcp)./abs(rc_nrw)) < 1e-3
% %             gam0_sign = 1;
% %         elseif max(abs(rc_nrw - rcm)./abs(rc_nrw)) < 1e-3
% %             gam0_sign = -1;
% %         else
%         nomatch = sum(gam0_sign==0);
%         if nomatch > 0
%             % If trans coefs don't match with either sign, warn and perform
%             % calculation with positive sign
%             error('Reverse reflection coefficient does not match NRW reflection coefficient')
% %             gam0_sign = 1;
%         end
%     else
%         gam0_sign = parser.Results.gam0_sign;
%     end
%     gam0 = gam0.*gam0_sign;

    rc = rev_ref_coef(nrw_out.mu,gam0,gam);
%     disp(tc(330:360))
    
    % calculate S-params
    out = nrw_out(:,1);
    deno = 1-(rc.^2.*tc.^2);
    out.s11 = R1.^2.*(rc.*(1-tc.^2)./deno);
    out.s21 = R1.*R2.*(tc.*(1-rc.^2)./deno);
    out.s12 = out.s21;
    out.s22 = R2.^2.*(rc.*(1-tc.^2)./deno);
    
    
    if parser.Results.modes==2
        lambda_c2 = parser.Results.lambda_c2;
        betas = parser.Results.betas;
        if lambda_c2==0 || isempty(betas)
            error('lambda_c2 and betas must be provided for higher-order modes')
        end
%         out.s11 = out.s11*(1+betas(1)
%   need to duplicate all calculations involving lambda_c
    elseif parser.Results.modes~=1
        error('Modes must be set to 1 or 2')
    end
    
    varargout = cell(1,nargout);
    varargout{1} = out;
    if nargout==3
        varargout{2} = gam_sign;
        varargout{3} = gam0_sign;
    elseif nargout==6
        varargout{2} = gam_sign;
        varargout{3} = gam0_sign;
        varargout{4} = tc_nrw;
        varargout{5} = tcp;
        varargout{6} = tcm;
    end
    
end