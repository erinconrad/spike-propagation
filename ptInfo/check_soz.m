function check_soz(pt)

for i = 1:31
    all_soz = [];
    if isempty(pt(i).newSOZChs) == 1, continue; end
    if ismember(i,[20,22,24,26]) == 1
        continue
    end
    
    for j = 1:length(pt(i).sz)
        for k = 1:length(pt(i).sz(j).electrodes)
            
            if isempty(pt(i).sz(j).chs) == 1, continue; end
            
            if pt(i).sz(j).chs(k) == 0, continue; end
            
            % Make sure it matches number
            if strcmp(pt(i).sz(j).electrodes{k},...
                    pt(i).electrodeData.electrodes(pt(i).sz(j).chs(k)).name) == 0
                error('What\n');
            end
            
            
            
        end
        
        % Add to array
        all_soz = [all_soz;pt(i).sz(j).chs];
    end
    
    all_soz = unique(all_soz);
    if isequal(all_soz,pt(i).newSOZChs) == 0
        error('What\n');
    end
    
end

end