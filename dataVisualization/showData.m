% This function plots the EEG voltages for a desired time and channels

function showData(Patient,pt,startTime,duration,chIds)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
dataName = Patient(pt).ieeg_name;
electrodeFile = Patient(pt).electrode_labels;
ptname = Patient(pt).name;

%% load seizure info and channel info
ptInfo = loadjson(jsonfile);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;
times = [startTime,startTime+duration];

[gdf,~,extraoutput] = getSpikeTimes(times,ptname,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,1,0,13,300);
values = extraoutput{1};
unignoredChLabels = extraoutput{2};
nchs = length(chIds);
whichCh = chIds(1:min(10,nchs));
plottimes =  [1:size(values,1)]/fs;

colors = {'b','r','g','c','m','b','r','g','c','m'};
figure
range = 0;
pl = zeros(length(whichCh),1);

for i = 1:length(whichCh)
    ch = whichCh(i);
    amps = values(:,ch) - range;
    pl(i) = plot(plottimes,amps,colors{i});
    hold on
    range = range + max(values(:,ch)) - min(values(:,ch));
end

legnames = unignoredChLabels(whichCh);
legend(pl,legnames,'Location','northeast');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 duration]);
title(sprintf('Data for %s time %d',...
    ptname,startTime));


end