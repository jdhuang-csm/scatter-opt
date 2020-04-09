function tc = nrw_trans_coef(s11,s21)
    import scatter_opt.*
	rc = nrw_ref_coef(s11,s21);
    tc = (s11+s21-rc)./(1-(s11+s21).*rc);
end