function pt = getSOZChs(pt)
    for i = 1:length(pt)
        for j = 1:length(pt(i).sz)
            chnames = pt(i).sz(j).electrodes;
            chnums = zeros(length(chnames),1);
            for ich = 1:length(chnames)
                for ich = 1:length(chnames)
                    [Lia,chIds] = ismember(chnames{ich},pt(i).electrodeData.unignoredChs);
                    if Lia == 0
                        fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                            chnames{ich},pt(i).name);
                        error('');
                        
                    end
                    chnums(ich) = chIds;
                end
                
                
                
            end
            pt(i).sz(j).chs = chnums;
        
        end
    end

end