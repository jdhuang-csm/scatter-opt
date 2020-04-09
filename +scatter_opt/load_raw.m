function raw = load_raw(file,varargin)
    parser = inputParser;
    addOptional(parser,'HeaderLines',0)
    addParameter(parser,'Delimiter',' ')
    addParameter(parser,'PhaseUnwrap',true)
    parse(parser,varargin{:})
    
    raw = readtable(file,'HeaderLines',parser.Results.HeaderLines,...
        'Delimiter',parser.Results.Delimiter);
    raw.Properties.VariableNames = {'freq','s11_db','s11_phase','s21_db','s21_phase','s12_db','s12_phase','s22_db','s22_phase'};
    
    s_ij = {'s11' 's21' 's12' 's22'};
    % phase unwrapping
    if parser.Results.PhaseUnwrap
        for i=1:4
            si = s_ij{i};
            phzdiff = diff(raw.(strcat(si,'_phase')));
            jumps = abs(phzdiff)>180;
            if sum(jumps)>0
                index = find(jumps);
                for idx = 1:length(index)
                    j = index(idx);
                    raw.(strcat(si,'_phase'))(j+1:end) = ...
                        raw.(strcat(si,'_phase'))(j+1:end) - sign(phzdiff(j))*360;
                end
            end
        end
    end
            
end