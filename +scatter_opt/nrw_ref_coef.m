function rc = nrw_ref_coef(s11,s21)
    X = (s11.^2 - s21.^2 + 1)./(2.*s11);
    rcp = X + sqrt(X.^2 - 1);
    rcm = X - sqrt(X.^2 - 1);
    % correct root given by |rc| < 1
    rc = zeros(length(s11),1);
    rc(abs(rcp)<1) = rcp(abs(rcp)<1);
    rc(abs(rcm)<1) = rcm(abs(rcm)<1);
    nomatch = find(rc==0,1);
    if ~isempty(nomatch)
        error('Neither sign gives |rc| < 1')
    end
%     if max(abs(rcp)) < 1
%         rc = rcp;
%     elseif max(abs(rcm)) < 1
%         rc = rcm;
%     end
end