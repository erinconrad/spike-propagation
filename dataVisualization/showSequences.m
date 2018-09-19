function showSequences(P,pts,whichSeq,nseq,ic)
% This is another function to plot sequences, using the spike times from
% the inputted structure

if isempty(whichSeq) == 0
    %ignore nseq
    nseq = length(whichSeq);
end

%% Parameters
prows = 2;
columns =  nseq/prows;
surroundtime = 2;

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);



for pt = pts

mkdir([resultsFolder,'plots/',P(pt).name,'/exampleSeqs/']);

if isfield(P(pt),'sz') == 0
    continue
end

if isfield(P(pt).sz(1),'seq_matrix') == 0
    continue
end

if isempty(P(pt).sz(1).seq_matrix) == 1
    continue
end


dataName = P(pt).ieeg_name;
ptname = P(pt).name;


%% Concatenate all sequences
allSeq = [];
for j = 1:length(P(pt).sz)
    allSeq = [allSeq,P(pt).sz(j).seq_matrix];
end

szTimes = [];
for j = 1:length(P(pt).sz)
    szTimes = [szTimes;P(pt).sz(j).onset P(pt).sz(j).offset];
end

firstSpikes = min(allSeq,[],1);
nonIctalSeq = allSeq(:,~any(firstSpikes >= szTimes(:,1) & ...
    firstSpikes <= szTimes(:,2),1));
ictalSeq = allSeq(:,any(firstSpikes >= szTimes(:,1) & ...
    firstSpikes <= szTimes(:,2),1));

%% If chunk not specified, randomly pick non-ictal sequences to plot
if isempty(whichSeq) == 1
    
    
    
    % Pick random set of non-ictal sequences
    if ic == 0
        fprintf('Picking a random set of non-ictal sequences from all seizures\n');
        y = randsample(size(nonIctalSeq,2),nseq);
        seqs = nonIctalSeq(:,y);

    
    elseif ic == 1
        fprintf('Picking a random set of ictal sequences from all seizures\n');
        y = randsample(size(ictalSeq,2),nseq);
        seqs = ictalSeq(:,y);
    end
else
    seqs = nonIctalSeq(:,whichSeq);
    
end



%% Define the structure that I will put plotting data into
seq = struct;

%% Loop through desired sequences and get the data
for i = 1:size(seqs,2)
    
    s = seqs(:,i);
    
    % Get the non nan times for the sequence
    spike_chs = find(isnan(s) == 0);
    spike_times = s(find(isnan(s) == 0));
    
    % sort them from earliest to latest times in the sequence
    [spike_times,I] = sort(spike_times);
    spike_chs = spike_chs(I);
    
    times = [spike_times(1)-surroundtime,spike_times(1)+surroundtime];
   
    %% Load EEG data info
    % calling this with 0 and 0 means I will just get basic info like sampling
    % rate and channel labels
    data = getiEEGData(dataName,0,0,pwfile);  
    fs = data.fs;

    %% Get the data for these times
    thresh.tmul = P(pt).thresh.tmul;
    thresh.absthresh = P(pt).thresh.absthresh;
    [gdf,extraoutput] = getSpikesSimple(P,pt,times,4,thresh);
    values = extraoutput.values;
    unignoredChLabels = P(pt).electrodeData.unignoredChs;
    plottimes =  [1:size(values,1)]/fs;

    %% Get the spike times
    
    spikes = [spike_chs,spike_times];
    
    %spikes = gdf;
    

    seq(i).seq = s;
    seq(i).spikes = spikes;
    seq(i).plottimes = plottimes;
    seq(i).values = values;
    seq(i).whichCh = spike_chs;
    seq(i).time_col = spike_times;
end

%% Plot 
colors = {'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m'};
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

for s = 1:nseq
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

    %{
    text(0.1,0.1,sprintf('Sequence %d %d s',seq(s).seq,...
        seq(s).time_col(1)),'units','normalized','FontSize',15);
    %}
    
    text(0.1,0.1,sprintf('Plot starts at %d s',seq(s).time_col(1)-surroundtime),...
        'units','normalized','FontSize',15);
    
    set(gca,'fontsize',15);
    xticks = 1:2:surroundtime*2-1;
    yl = ylim;
    yloc = yl(1)+(yl(2)-yl(1))*0.05;
    for t = 1:length(xticks)
       tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
    end

end



if isempty(whichSeq) == 0
    outputFile = [ptname,'_sequences'];
elseif isempty(whichSeq) == 1
    if ic == 1
        outputFile = [ptname,'_sequences_ic'];
    elseif ic == 0
        outputFile = [ptname,'_sequences_inter'];
    end
end

saveas(gcf,[resultsFolder,'plots/',P(pt).name,'/exampleSeqs/',outputFile,'.png'])


%% Plot again but make it pretty
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

for s = 1:nseq
    whichcol = ceil(s/2);
    whichrow = mod(s+1,2)+1;
    subplot('Position',[(whichcol-1)*1/columns (prows-whichrow)*1/prows 1/columns 1/prows])
    range = 0;
    pl = zeros(length(seq(s).whichCh),1);
    for i = 1:length(seq(s).whichCh)
        ch = seq(s).whichCh(i);
        amps = seq(s).values(:,ch) - range;
        pl(i) = plot(seq(s).plottimes,amps,'k');
        hold on


        % find the times of the spikes with the desired channel
        spiketimes = seq(s).spikes(seq(s).spikes(:,1) == ch,2)-seq(s).time_col(1)+surroundtime;

        spikeamp = ones(size(spiketimes,1),1)*max(seq(s).values(:,ch))-range;

        scatter(spiketimes,spikeamp,80,'k','filled');
        range = range + max(seq(s).values(:,ch)) - min(seq(s).values(:,ch));
    end
    legnames = unignoredChLabels(seq(s).whichCh);
    %legend(pl,legnames,'Location','northeast');
    xlabel('Time (s)');
    %ylabel('Amplitude');
    set(gca,'YTickLabel',[]);

    %{
    text(0.1,0.1,sprintf('Sequence %d %d s',seq(s).seq,...
        seq(s).time_col(1)),'units','normalized','FontSize',15);
    %}
    
    %text(0.1,0.1,sprintf('Plot starts at %d s',seq(s).time_col(1)-surroundtime),...
    %    'units','normalized','FontSize',15);
    
    set(gca,'fontsize',15);
    xticks = 1:2:surroundtime*2-1;
    yl = ylim;
    yloc = yl(1)+(yl(2)-yl(1))*0.05;
    %{
    for t = 1:length(xticks)
       tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
    end
    %}

end




outputFile = [outputFile,'_pretty'];

saveas(gcf,[resultsFolder,'plots/',P(pt).name,'/exampleSeqs/',outputFile,'.png'])

end


end