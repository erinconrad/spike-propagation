function visualizeSequences(Patient,whichSz,block,s)


%% Parameters to change every time

doplot =  1;

% data name (for ieeg.org)
dataName = 'HUP80_phaseII';
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
csvFile = 'HUP080_T1_19991213_electrode_labels.csv';
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
ptname = 'HUP080';
%ptname = 'HUP078';

% The number of the patient
pt = 80;
%pt = 78;

% Remove EKG artifact?
rmEKGArtifact = 0;


%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);

Patient(pt).seizures = ptInfo.PATIENTS.(ptname).Events.Ictal;
szNames = fieldnames(Patient(pt).seizures);

%% Define seizure onset and offset times for each seizure
for i = 1:length(fieldnames(Patient(pt).seizures))
    Patient(pt).sz(i).onset = Patient(pt).seizures.(szNames{i}).SeizureEEC;
    Patient(pt).sz(i).offset = Patient(pt).seizures.(szNames{i}).SeizureEnd;
    Patient(pt).sz(i).electrodes = {};
    for j = 1:length(Patient(pt).seizures.(szNames{i}).SEIZURE_ONSET_ELECTRODES)
        Patient(pt).sz(i).electrodes = [Patient(pt).sz(i).electrodes;Patient(pt).seizures.(szNames{i}).SEIZURE_ONSET_ELECTRODES{j}];
        
    end
    
end

%% Define the start and stop times of each seizure
blockOnset = Patient(80).sz(whichSz).runTimes(block,1);
time_col = Patient(pt).sz(whichSz).block(block).data.sequences(:,(s-1)*2+2);
chan_col = Patient(pt).sz(whichSz).block(block).data.sequences(:,(s-1)*2+1);
times = blockOnset+[time_col(1)-8,time_col(1)+8];
whichCh = chan_col(1:5);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% calculate gdf (spike times and locations) and output the data in that time
fprintf('Detecting spikes\n');
[gdf,~,extraoutput] = getSpikeTimes(times,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0);
values = extraoutput{1};
unignoredChLabels = extraoutput{2};
plottimes =  [1:size(values,1)]/fs;


%% calculate gdf and values of EKG channels
[gdfEKG,~,extraoutputEKG] = getSpikeTimes(times,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,1);
valuesEKG = extraoutputEKG{1};
unignoredChLabelsEKG = extraoutputEKG{2};
plottimesEKG =  [1:size(values,1)]/fs;

%% remove spikes that occur too close to EKG channel spikes
if rmEKGArtifact == 1
    prox = 0.01; %10 ms
    oldgdf =  gdf;
    gdf = removeEKGArtifact(gdf,gdfEKG,prox);
end


%% Plot 
if doplot == 1
colors = {'b','r','g','c','m'};
figure
range = 0;
pl = zeros(length(whichCh),1);
for i = 1:length(whichCh)
    ch = whichCh(i);
    amps = values(:,ch) - range;
    pl(i) = plot(plottimes,amps,colors{i});
    hold on
    
    
    % find the times of the spikes with the desired channel
    spiketimes = gdf(gdf(:,1) == ch,2);
    
    spikeamp = ones(size(spiketimes,1),1)*max(values(:,ch))-range;
    
    scatter(spiketimes,spikeamp,colors{i});
    range = range + max(values(:,ch)) - min(values(:,ch));
end
legnames = unignoredChLabels(whichCh);
legend(pl,legnames,'Location','northeast');
xlabel('Time (s)');
ylabel('Amplitude');

title(sprintf('Data and spike detections for %s seizure %d block %d sequences %d',...
    ptname,whichSz,block,s));
set(gca,'fontsize',15);




if 1 == 0
% add on EKG channel
yl = ylim;
ampsEKG = valuesEKG(:,1)-range;
plot(plottimesEKG,ampsEKG,'k');
spiketimesEKG =  gdfEKG(gdfEKG(:,1) == 1,2);
spikeampEKG = ones(size(spiketimesEKG,1),1)*max(valuesEKG(:,1))-range;
scatter(spiketimesEKG,spikeampEKG,'k');
yl(1) =  yl(1)-range*2;

ylim(yl)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.6, 1, 0.4]);

outputFile = [ptname,'_','_sz_',sprintf('%d',whichSz),'_block_',...
    sprintf('%d',block),'seq_',sprintf('%d',s),'.png'];

saveas(gcf,[resultsFolder,outputFile])


fprintf('no\n');

end

end