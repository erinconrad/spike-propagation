dataName = 'I022_P001_D05';

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
data = getiEEGData(dataName,0,0,pwfile);

vetted = load('/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/realSpikes/vettedvar.mat');
firstSpike = vetted.vetted{1}.annots(1);
%startAndEnd = round([firstSpike.start-100*1e6,firstSpike.stop+100*1e6]/1e6*data.fs);


%visualizeSpikes(dataName,[],'Cog patient',80,times,whichCh)