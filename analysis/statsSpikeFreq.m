function statsSpikeFreq(pt,whichPt,window,nwindows)


%% Remove EKG artifact and depth electrodes
rmEKG = 1;
rmDepth = 1;


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
gdfFolder = [resultsFolder,'gdf/'];
nchs = length(pt(whichPt).channels);
gdf_all{1} = [];

for j = 1:length(pt(whichPt).sz)
    
    
end

end