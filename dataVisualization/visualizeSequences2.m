function visualizeSequences2(Patient,pt,whichSz,block,whichSeq)
% This is another function to plot sequences, using the spike times from
% the inputted structure

%% Parameters to change every time
prows = 2;
columns =  length(whichSeq)/prows;

surroundtime = 4;

dummyRun = 0;
vanleer = 0;
vtime = 0;
outputData = 1;
keepEKG = 0;
ignore = 1;
funnyname = 0;
onlyclean = 0;
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
%pt = 78;

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


for i = 1:length(whichSeq)
    
    s = whichSeq(i);
    if isnan(block) == 1
        
        onsettime = Patient(pt).sz(whichSz).runTimes(1);
        
        time_col = Patient(pt).sz(whichSz).data.sequences(:,(s-1)*2+2);
        chan_col = Patient(pt).sz(whichSz).data.sequences(:,(s-1)*2+1);
        times = onsettime + [time_col(1)-surroundtime,time_col(1)+surroundtime];
        whichCh = chan_col(1:5);
    else
    
        blockOnset = Patient(pt).sz(whichSz).runTimes(block,1);
    
        if onlyclean == 0
            time_col = Patient(pt).sz(whichSz).block(block).data.sequences(:,(s-1)*2+2);
            chan_col = Patient(pt).sz(whichSz).block(block).data.sequences(:,(s-1)*2+1);
        elseif onlyclean == 1
            time_col = Patient(pt).sz(whichSz).block(block).data.cleanseq(:,(s-1)*2+2);
            chan_col = Patient(pt).sz(whichSz).block(block).data.cleanseq(:,(s-1)*2+1);
        end
        times = blockOnset+[time_col(1)-surroundtime,time_col(1)+surroundtime];
        whichCh = chan_col(1:5);
    end

    %% Load EEG data info
    % calling this with 0 and 0 means I will just get basic info like sampling
    % rate and channel labels
    data = getiEEGData(dataName,0,0,pwfile);  
    fs = data.fs;

    %% Get the data for these times
    fprintf('Detecting spikes\n');
    [~,~,extraoutput] = getSpikeTimes(times,dataName,electrodeFile,ptInfo,pwfile,...
        dummyRun,vanleer,vtime,outputData,keepEKG,ignore,funnyname);
    values = extraoutput{1};
    unignoredChLabels = extraoutput{2};
    plottimes =  [1:size(values,1)]/fs;

    %% Get the spike times
    spikes = [chan_col,time_col];
    
    seq(i).seq = s;
    seq(i).spikes = spikes;
    seq(i).plottimes = plottimes;
    seq(i).values = values;
    seq(i).whichCh = whichCh;
    seq(i).time_col = time_col;

end

%% Plot 
colors = {'b','r','g','c','m'};
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.9, 0.8]);

for s = 1:length(whichSeq)
    whichcol = ceil(s/2);
    whichrow = mod(s+1,2)+1;
    subplot('Position',[(whichcol-1)*1/columns (prows-whichrow)*1/prows 1/columns 1/prows])
    range = 0;
    pl = zeros(length(seq(s).whichCh),1);
    for i = 1:length(seq(s).whichCh)
        ch = seq(s).whichCh(i);
        amps = seq(s).values(:,ch) - range;
        pl(i) = plot(seq(s).plottimes,amps,colors{i});
        hold on


        % find the times of the spikes with the desired channel
        spiketimes = seq(s).spikes(seq(s).spikes(:,1) == ch,2)-seq(s).time_col(1)+surroundtime;

        spikeamp = ones(size(spiketimes,1),1)*max(seq(s).values(:,ch))-range;

        scatter(spiketimes,spikeamp,60,colors{i},'filled');
        range = range + max(seq(s).values(:,ch)) - min(seq(s).values(:,ch));
    end
    legnames = unignoredChLabels(seq(s).whichCh);
    legend(pl,legnames,'Location','northeast');
    xlabel('Time (s)');
    %ylabel('Amplitude');
    set(gca,'YTickLabel',[]);

    text(0.7,0.1,sprintf('Sequence %d',seq(s).seq),'units','normalized','FontSize',15);
    set(gca,'fontsize',15);

end


%{
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
%}


outputFile = [ptname,'_','_sz_',sprintf('%d',whichSz),'_block_',...
    sprintf('%d',block),'.png'];

saveas(gcf,[resultsFolder,outputFile])


fprintf('no\n');


end