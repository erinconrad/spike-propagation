function pt = newSzOnsetChs(pt)

for i = 1:length(pt)
    szElecs = [];
    for j = 1:length(pt(i).sz)

        for k = 1:length(pt(i).sz(j).electrodes)
                szElecs = [szElecs,pt(i).sz(j).electrodes{k}];
        end
    end
    
    szElecs = unique(szElecs);
    pt(i).newSOZNames = szElecs;
    
    
    chnums = zeros(length(szElecs),1);
    for ich = 1:length(szElecs)
        for ich = 1:length(szElecs)
            [Lia,chIds] = ismember(szElecs{ich},pt(i).electrodeData.unignoredChs);
            if Lia == 0
                fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                    szElecs{ich},pt(i).name);
                error('look');

            end
            chnums(ich) = chIds;
        end



    end
    
    pt(i).newSOZChs = chnums;
    
end



end