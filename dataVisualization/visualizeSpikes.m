function visualizeSpikes(dataName,csvFile,ptname,pt,times,whichCh)


%% Parameters to change every time

% ictal or pre-ictal
ictal = 1;
doplot =  1;
% which seizure
whichSz = 2;

if nargin == 0
% which channels
%whichCh = [50,51,52];





% data name (for ieeg.org)
dataName = 'I022_P001_D05';
%dataName = 'HUP80_phaseII';
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
csvFile = [];
%csvFile = 'HUP080_T1_19991213_electrode_labels.csv';
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
ptname = 'Cog patient';
%ptname = 'HUP080';
%ptname = 'HUP078';

% The number of the patient
pt = 80;
%pt = 78;

times = [1616,1636];
end

% Remove EKG artifact?
rmEKGArtifact = 0;


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
electrodeFile = 0;
ptInfo = 1;
if isempty(csvFile) == 0
%% Get paths and load seizure info and channel info

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

duration = Patient(pt).sz(whichSz).offset - Patient(pt).sz(whichSz).onset;

if ictal == 1
    times = [Patient(pt).sz(whichSz).onset,Patient(pt).sz(whichSz).offset];
else
    times = [Patient(pt).sz(whichSz).onset-duration*1,Patient(pt).sz(whichSz).offset];
   
end
end

%times = 34117+[800,850];
%times = 36000 + [0,500];
%times = [37000 37116];


%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% calculate gdf (spike times and locations) and output the data in that time
fprintf('Detecting spikes\n');
[gdf,~,extraoutput] = getSpikeTimes(times,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,0,1);
values = extraoutput{1};
unignoredChLabels = extraoutput{2};
plottimes =  [1:size(values,1)]/fs;




%% remove spikes that occur too close to EKG channel spikes
if rmEKGArtifact == 1
    %% calculate gdf and values of EKG channels
[gdfEKG,~,extraoutputEKG] = getSpikeTimes(times,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,1,0,1);
valuesEKG = extraoutputEKG{1};
unignoredChLabelsEKG = extraoutputEKG{2};
plottimesEKG =  [1:size(values,1)]/fs;
    
    
    prox = 0.01; %10 ms
    oldgdf =  gdf;
    gdf = removeEKGArtifact(gdf,gdfEKG,prox);
end

%% Get which channels to plot
if isempty(csvFile) == 0
[~,chIds] = ismember(Patient(pt).sz(whichSz).electrodes,unignoredChLabels);
whichCh = chIds(1:5);
end

%% Plot 
if doplot == 1
colors = {'b','r','g','c','m','b','r','g','c','m'};
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

if ictal == 1, ictext = 'Ictal'; else, ictext = 'Pre-ictal and ictal'; end
title(sprintf('%s data and spike detections for %s seizure %d',...
    ictext,ptname,whichSz));
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

outputFile = [ptname,'_',ictext,'_sz_',sprintf('%d',whichSz),'.png'];

%saveas(gcf,[resultsFolder,outputFile])


fprintf('no\n');

end

end