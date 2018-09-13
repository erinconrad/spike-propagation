function portVisualizeSequences(P,pts,whichSz,whichSeq,allSpikes)
% This is another function to plot sequences, using the spike times from
% the inputted structure

for pt = pts

%% Parameters to change every time
prows = 2;
columns =  6/prows;

surroundtime = 7.5;

dummyRun = 0;
vanleer = 0;
vtime = 0;
outputData = 1;
keepEKG = 0;
ignore = 1;
funnyname = 0;
onlyclean = 0;
% data name (for ieeg.org)
dataName = P(pt).ieeg_name;
%dataName = 'HUP78_phaseII-Annotations';  

% CSV file with electrode locations
%csvFile = 'HUP078_T1_19971218_electrode_labels.csv';

% The patient name with format as used in the json file
ptname = P(pt).name;
%ptname = 'HUP078';

% The number of the patient
%pt = 78;

%% Find the ictal chunks so that I can ignore them
ictal = zeros(size(P(pt).sz(whichSz).timesRL,1),1);
sztimes = [P(pt).sz(whichSz).onset, P(pt).sz(whichSz).offset];
for k = 1:size(P(pt).sz(whichSz).timesRL,1)
    chunkTimes = P(pt).sz(whichSz).timesRL + P(pt).sz(whichSz).runTimes(1,1);
    
    if chunkTimes(k,1) < sztimes(2) && ...
            chunkTimes(k,2) > sztimes(1)
        ictal(k) = 1;
        
    end
    
end

%% If chunk not specified, pseudorandomly pick 6 non-ictal sequences to plot
if isempty(whichSeq) == 1
    % first pick 6 random chunks that don't have ictal data
    allchunks = 1:size(P(pt).sz(whichSz).timesRL,1);
    nonictal = allchunks(~ictal);
    
    fprintf('in loop\n');
    while 1
        specchunks = randsample(nonictal,6);
        seqok = zeros(6,1);
        for i = 1:length(specchunks)
            thing = P(pt).sz(whichSz).blockRL(specchunks(i));
            if isfield(thing,'sIdx') == 1
                if isempty(thing.sIdx) == 0
                    seqok(i) = 1;
                end
            end
        end
        if sum(seqok) == 6
            break
        end
        
    end
    
    fprintf('done with loop\n');
    % pick a sequence from each chunk to plot
    whichSeq = zeros(6,1);
    for i = 1:length(specchunks)
         whichSeq(i) = randsample(P(pt).sz(whichSz).blockRL(specchunks(i)).sIdx,1);
    end
    
end

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);


if allSpikes == 0, asText = 'justSeq'; elseif allSpikes ==1, asText = 'newSpike'; end

%% Define the start and stop times of each seizure

for i = 1:length(whichSeq)
    
    s = whichSeq(i);
   
        
    onsettime = P(pt).sz(whichSz).runTimes(1);

    time_col = P(pt).sz(whichSz).data.sequences(:,(s-1)*2+2);
    chan_col = P(pt).sz(whichSz).data.sequences(:,(s-1)*2+1);
    times = [time_col(1)-surroundtime,time_col(1)+surroundtime];
    %times = onsettime + [time_col(1)-surroundtime,time_col(1)+surroundtime];
    whichCh = chan_col(1:5);
   
    %% Load EEG data info
    % calling this with 0 and 0 means I will just get basic info like sampling
    % rate and channel labels
    data = getiEEGData(dataName,0,0,pwfile);  
    fs = data.fs;

    %% Get the data for these times
    fprintf('Detecting spikes\n');
    thresh.tmul = P(pt).thresh.tmul;
    thresh.absthresh = P(pt).thresh.absthresh;
    [gdf,extraoutput] = getSpikesSimple(P,pt,times,4,thresh);
    %[gdf,~,extraoutput] = getSpikeTimes(times,ptname,dataName,1,ptInfo,pwfile,...
    %    dummyRun,vanleer,vtime,outputData,keepEKG,ignore,funnyname,8,300);%P(pt).tmul,P(pt).absthresh);
    values = extraoutput.values;
    unignoredChLabels = P(pt).electrodeData.unignoredChs;
    plottimes =  [1:size(values,1)]/fs;

    %% Get the spike times
    if allSpikes == 0
        spikes = [chan_col,time_col];
    elseif allSpikes == 1
        spikes = gdf;
        spikes(:,2) = spikes(:,2) + time_col(1) - surroundtime;
    end

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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

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

        scatter(spiketimes,spikeamp,80,colors{i},'filled');
        range = range + max(seq(s).values(:,ch)) - min(seq(s).values(:,ch));
    end
    legnames = unignoredChLabels(seq(s).whichCh);
    legend(pl,legnames,'Location','northeast');
    xlabel('Time (s)');
    %ylabel('Amplitude');
    set(gca,'YTickLabel',[]);

    text(0.1,0.1,sprintf('Sequence %d %d s',seq(s).seq,...
        seq(s).time_col(1)),'units','normalized','FontSize',15);
    
    text(0.1,0.9,sprintf('Plot starts at %d s',seq(s).time_col(1)-surroundtime),...
        'units','normalized','FontSize',15);
    
    set(gca,'fontsize',15);
    xticks = 1:2:surroundtime*2-1;
    yl = ylim;
    yloc = yl(1)+(yl(2)-yl(1))*0.05;
    for t = 1:length(xticks)
       tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
    end

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


outputFile = [ptname,'_sz_',sprintf('%d',whichSz),'_',asText,'_tmul8_.png'];

saveas(gcf,[resultsFolder,'plots/',P(pt).name,'/',outputFile])


end


end