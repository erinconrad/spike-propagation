function addThresh(whichPts)

[~,~,scriptFolder,resultsFolder,~] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];
gdfFile = 'long_gdf.mat';
load(structFolder,gdfFile);

for whichPts = whichPt
    name = pt(whichPt).name;
    [~,~,thresh,dmin,~] = ieegAndElectodeNames(name);
    pt(whichPt).thresh = thresh;
    pt(whichPt).dmin = dmin;
    
end

save(structFolder,gdfFile);


end