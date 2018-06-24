% This function plots the EEG voltages for a desired time and channels

function showData(Patient,pt,startTime,duration,chIds)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
dataName = Patient(pt).ieeg_name;
electrodeFile = Patient(pt).electrode_labels;
ptname = Patient(pt).name;
fs = Patient(pt).fs;

%% load seizure info and channel info
ptInfo = loadjson(jsonfile);

%% Load EEG data info
times = [startTime,startTime+duration];
[gdf,extraOutput] = getSpikesSimple(Patient,pt,times,4);
values = extraOutput.values;
unignoredChLabels = Patient(pt).electrodeData.unignoredChs;


%{
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;
times = [startTime,startTime+duration];

[gdf,~,extraoutput] = getSpikeTimes(times,ptname,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,1,0,13,300,4);
values = extraoutput{1};
unignoredChLabels = extraoutput{2};
%}
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
    
    nbins = ceil(length(amps)/60/512);
    noise_bin = zeros(nbins,1);
    col = zeros(nbins,3);
    for j = 1:nbins
       noise_bin(j) = findNoisyPeriods(values((j-1)*60*512+1:...
           min(length(amps),j*60*512),ch),2); 
       if noise_bin(j) == 1, col(j,:) = [1 0 0]; else, col(j,:) = [0 0 1]; end 
    end
    old_range = range;
    range = range + max(values(:,ch)) - min(values(:,ch));
    noise_bin = noise_bin*(max(values(:,ch)) - min(values(:,ch)))*3/4 - old_range*2;
    
    
    
    scatter(linspace(plottimes(1),plottimes(end),nbins),noise_bin*0.5,100,col,'filled');
    
    
    
end

legnames = unignoredChLabels(whichCh);
legend(pl,legnames,'Location','northeast');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0 duration]);
title(sprintf('Data for %s time %d',...
    ptname,startTime));


end