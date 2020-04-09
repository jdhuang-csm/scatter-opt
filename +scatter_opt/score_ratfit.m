function score = score_ratfit(fit,actual,freq,varargin)
    parser = inputParser;
    addOptional(parser,'Part','both')
    parse(parser,varargin{:})
    
    [pred,outfreq] = freqresp(fit,freq);
    
    if strcmp(parser.Results.Part,'real')
        pred = real(pred);
        actual = real(actual);
    elseif strcmp(parser.Results.Part,'imag')
        pred = imag(pred);
        actual = imag(actual);
    elseif ~strcmp(parser.Results.Part,'both')
        error('Invalid Part argument. Options are both, imag, real')
    end
    
    resid = pred - actual;
    %rwr = real(resid)./real(pred);
    %iwr = imag(resid)./imag(pred);
    %disp(rwr.^2 + iwr.^2)
    % residual sum of squares
    rss = sum(real(resid).^2 + imag(resid).^2);
    %ssreg = sum((real(pred)-mean(real(actual))).^2) + sum((imag(pred)-mean(imag(actual))).^2);
    %rss = sum(rwr.^2 + iwr.^2);
    % total sum of squares
    ssr = sum((real(actual)-mean(real(actual))).^2);
    ssi = sum((imag(actual)-mean(imag(actual))).^2);
    sstot = ssr + ssi;
    score = 1 - (rss/sstot);
end