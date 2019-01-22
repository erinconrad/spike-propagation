function pt = updateSOZ(pt)

%HUP116
pt(20).sz(1).electrodes = {'RA1','RA2','RH1','RH2','RH3','RH4'};
pt(20).sz(2).electrodes = {'RA1','RA2','RH1','RH2','RH3'};
pt(20).sz(3).electrodes = {'RCC2'};

%Study016
pt(22).sz(1).electrodes = {'ROF4','RTG9','RTG21','RTG22','RTG23'};
pt(22).sz(2).electrodes = {'ROF2','ROF3','ROF4','RFG8','RFG9','RFG15'};
pt(22).sz(3).electrodes = {'ROF4','RFG18','RFG24'};
pt(22).sz(4).electrodes = {'ROF1','ROF2','ROF3','ROF4'};
pt(22).sz(5).electrodes = {'ROF1','ROF2','ROF3','ROF4'};
pt(22).sz(6).electrodes = {'ROF2','ROF3','ROF4','RFG9','RFG9'};

%Study019
pt(24).sz(1).electrodes = {'LT20','LT21','PT3','PT4','LT15','LT16'};
pt(24).sz(2).electrodes = {'PT3','PT4','LT15','LT16'};
pt(24).sz(3).electrodes = {'PT3','PT4','LT15','LT16','LT21','LT22','LT27',...
    'LT28','LT32','LT33'};
pt(24).sz(4).electrodes = {'PT3','PT4','LT15','LT16','LT21','LT22','LT27',...
    'LT28','LT32','LT33'};
pt(24).sz(5).electrodes = {'PT3','PT4','LT15','LT16','LT21','LT22','LT27',...
    'LT28','LT32','LT33'};

for i = [20 22 24]
    szElecs = {};
    for j = 1:length(pt(i).sz)

        for k = 1:length(pt(i).sz(j).electrodes)
            szElecs = [szElecs,pt(i).sz(j).electrodes{k}];
        end
    end
    
    szElecs = unique(szElecs);
    pt(i).newSOZNames = szElecs;
    
    
    chnums = zeros(length(szElecs),1);
    for ich = 1:length(szElecs)
        [Lia,chIds] = ismember(szElecs(ich),pt(i).electrodeData.unignoredChs);
        if Lia == 0
            fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                szElecs(ich),pt(i).name);
            error('look');

        end
        chnums(ich) = chIds;
    end
    
    pt(i).newSOZChs = chnums;
    
end

end

