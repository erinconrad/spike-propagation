function [electrodeFile,jsonfile,scriptFolder,resultsFolder] = fileLocations 


mainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/';
electrodeFile = [mainFolder,'data/','electrodeData/','HUP078_T1_19971218_electrode_labels.csv'];
jsonfile = [mainFolder,'data/','patientData/','DATA.json'];
scriptFolder = [mainFolder,'scripts/','my scripts/'];
resultsFolder = [mainFolder,'results/'];

end