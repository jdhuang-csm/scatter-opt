function conv = loadconvert_raw(file,varargin)
% load raw s_ij data and convert to complex
    import scatter_opt.*
    parser = inputParser;
    addOptional(parser,'HeaderLines',0)
    addOptional(parser,'Delimiter',' ')
    parse(parser,varargin{:})
    
    raw = load_raw(file,'HeaderLines',parser.Results.HeaderLines,...
        'Delimiter',parser.Results.Delimiter);
    s_cols = {'s11' 's21' 's12' 's22'};
    conv = raw(:,1);
    for i=1:4
        si = s_cols{i};
        mag = 10.^(raw.(strcat(si,'_db'))./20);
        real = mag.*cosd(raw.(strcat(si,'_phase')));
        imag = mag.*sind(raw.(strcat(si,'_phase')));
        conv.(si) = real + 1i*imag;
    end 
   
end