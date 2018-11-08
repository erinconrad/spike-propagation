function pt = getSOZChs(pt,whichChs)

for i = whichChs

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

for k = 1:length(ptnames)
    if strcmp(pt(i).name,ptnames{k}) == 1
        info = ptInfo.PATIENTS.(ptnames{k});
        szs = fieldnames(info.Events.Ictal);
        for j = 1:length(szs)
            sz = info.Events.Ictal.(szs{j});
            for whichSz = 1:length(pt(i).sz)
               if pt(i).sz(whichSz).onset == sz.SeizureEEC 
                    pt(i).sz(whichSz).electrodes = sz.SEIZURE_ONSET_ELECTRODES;
                   
               end
            end
            
        end
        
    end
    
end


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