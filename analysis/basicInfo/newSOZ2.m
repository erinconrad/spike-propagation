function pt = newSOZ2(pt,whichPts)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

for whichPt = whichPts
    name = pt(whichPt).name;
    pt(whichPt).newSOZChs = [];
    pt(whichPt).newSOZNames ={};
    pt(whichPt).newSzTimes = [];
    for i = 1:length(ptnames)
       if strcmp(name,ptnames{i}) == 1 ||  (strcmp(name,'HUP111') == 1 &&...
                   strcmp(ptnames{i},'HUP111A') ==1 )
            info = ptInfo.PATIENTS.(ptnames{i});
            pt(whichPt).outcome = info.Outcome;
            pt(whichPt).SeizureOnset = info.SeizureOnset;
            pt(whichPt).AgeOnset = info.AgeOnset;
            szs = fieldnames(info.Events.Ictal);
            for j = 1:length(szs)
                sz = info.Events.Ictal.(szs{j});
                pt(whichPt).newSzTimes = [pt(whichPt).newSzTimes;...
                    sz.SeizureEEC sz.SeizureEnd];
                for k = 1:length(sz.SEIZURE_ONSET_ELECTRODES)
                    pt(whichPt).newSOZNames = [pt(whichPt).newSOZNames,...
                    sz.SEIZURE_ONSET_ELECTRODES{k}];
                end
                
            end
           
       end
        
    end
    
    pt(whichPt).newSOZNames = unique(pt(whichPt).newSOZNames);
    szElecs = pt(whichPt).newSOZNames;
    chnums = zeros(length(szElecs),1);
    for ich = 1:length(szElecs)
        [Lia,chIds] = ismember(szElecs(ich),pt(whichPt).electrodeData.unignoredChs);
        if Lia == 0
            fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                szElecs(ich),pt(whichPt).name);
            %error('look');

        end
        chnums(ich) = chIds;
    end
    
    pt(whichPt).newSOZChs = chnums;
    
    
end




end