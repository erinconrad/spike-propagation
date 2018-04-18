function portVisualizeSpikes(Patient,pt,whichSz,ictal,szOnsetZone,chIds,tmul,absthresh)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

if ictal == 1, ictext = 'Ictal'; elseif ictal == 2, ictext = 'Pre-ictal';...
elseif ictal == 3, ictext = 'inter-ictal'; end
outputFolder = [resultsFolder,'spike verification/',Patient(pt).name,'/'];

if szOnsetZone == 0
    whichCh = chIds;
end

dataName = Patient(pt).ieeg_name;
electrodeFile = Patient(pt).electrode_labels;
ptname = Patient(pt).name;



%% load seizure info and channel info
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
if ictal == 1
    times = [Patient(pt).sz(whichSz).onset,Patient(pt).sz(whichSz).onset + 15];
elseif ictal == 2
    times = [Patient(pt).sz(whichSz).onset-15,Patient(pt).sz(whichSz).onset];
elseif ictal == 3
    times = [Patient(pt).sz(whichSz).onset-60*10,...
        Patient(pt).sz(whichSz).onset-60*10+15];
   
end

if times(1) < 0
    fprintf('funny seizure times for patient %s\n',ptname);
    return
end

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% calculate gdf (spike times and locations) and output the data in that time
fprintf('Detecting spikes\n');
[gdf,~,extraoutput] = getSpikeTimes(times,ptname,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,1,0,tmul,absthresh);
values = extraoutput{1};
unignoredChLabels = extraoutput{2};
tmul = extraoutput{3};
plottimes =  [1:size(values,1)]/fs;

if szOnsetZone == 1, elec_text = 'seizure onset zone'; else, elec_text = 'random'; end
if szOnsetZone == 1, elecf_text = 'SOZChs'; else, elecf_text = 'randomChs'; end


if isempty(gdf) == 1
    fprintf('No spikes detected for patient %s, time %s, with threshold %s\n',...
        ptname,ictext,tmul);
    figure
    plot(1,1)
    
    title(sprintf('%s data and spike detections for %s seizure %d\n showing %s electrodes, threshold of %d',...
    ictext,ptname,whichSz,elec_text,tmul));

    outputFile = [ptname,'_',ictext,'_sz_',sprintf('%d',whichSz),elecf_text,'_threshold_',sprintf('%d',tmul),'.png'];

    saveas(gcf,[outputFolder,outputFile])
    close(gcf);
    return
end

%% Get which channels to plot
if szOnsetZone == 1
    [~,chIds] = ismember(Patient(pt).sz(whichSz).electrodes,unignoredChLabels);
    nchs = length(chIds);
    whichCh = chIds(1:min(10,nchs));
end

%% Plot 

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

title(sprintf('%s data and spike detections for %s seizure %d\n showing %s electrodes, threshold of %d, absthresh of %d',...
    ictext,ptname,whichSz,elec_text,tmul,absthresh));
set(gca,'fontsize',15);



set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.6, 1, 0.6]);


outputFile = [ptname,'_',ictext,'_sz_',sprintf('%d',whichSz),elecf_text,'_threshold_',sprintf('%d',tmul),'_absthresh_',sprintf('%d',absthresh),'.png'];

saveas(gcf,[outputFolder,outputFile])
close(gcf);

fprintf('\n');



end