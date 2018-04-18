function P = addElectrodeData(P)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
ptInfo = loadjson(jsonfile);
filename = 'pWithElectrodes.mat';

for i = 1:length(P)
    
    name = P(i).name;
    [ieeg_name,electrode_name] =  ieegAndElectodeNames(name);
    electrodeFile = [electrodeFolder,electrode_name];
    
    if isempty(ieeg_name) == 1 || isempty(electrode_name) == 1
       continue
    end
    
    dummyRun = 1;
    outputData = 0;
    [~,electrodeData,~] = getSpikeTimes(0,name,ieeg_name,electrodeFile,ptInfo,pwfile,...
    dummyRun,0,0,outputData,0,1,0);
    P(i).electrodeData = electrodeData;
    
end

save([resultsFolder,filename],'P');

end