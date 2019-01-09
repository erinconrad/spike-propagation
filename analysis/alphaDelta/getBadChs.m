function [badChNums,badChNamesOut] = getBadChs(pt,whichPt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

for i = 1:length(ptnames)
    name = ptnames{i};
    info = ptInfo.PATIENTS.(ptnames{i});
    if strcmp(name,pt(whichPt).name) == 1
        badChs = info.IGNORE_ELECTRODES;
        
    end
    
end

badChNums = [];
badChNamesOut = {};

for i = 1:length(badChs)
    for j = 1:length(pt(whichPt).electrodeData.electrodes)
        if strcmp(pt(whichPt).electrodeData.electrodes(j).name,badChs{i}) == 1
            badChNums = [badChNums;j];
            badChNamesOut = [badChNamesOut; pt(whichPt).electrodeData.electrodes(j).name];
        end
    end
end

end