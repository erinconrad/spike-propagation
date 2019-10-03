function showSequences(P,pts,whichSeq,nseq,ic,out_folder,which_cluster,example)
% This is another function to plot sequences, using the spike times from
% the inputted structure

save_plots = 0;

if isempty(pts) == 1
    for i = 1:length(P)
        if isempty(P(i).seq_matrix) == 0
            pts = [pts,i];
        end
    end
end

if isempty(whichSeq) == 0
    %ignore nseq
    nseq = size(whichSeq,2);
end

%% Figure out how many figures I need
nfigs = ceil(nseq/10);
seq_fig = cell(nfigs,1);
for i = 1:length(seq_fig)
    seq_fig{i} = (i-1)*10+1:min(i*10,nseq);
end
prows = 2;
columns =  5;
surroundtime = 2;

%% Get paths and load seizure info and channel info
if example == 0
    [electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
    p1 = genpath(scriptFolder);
    addpath(p1);
    ptInfo = loadjson(jsonfile);
end



for pt = pts

if isempty(out_folder) ==1
    out_folder = [resultsFolder,'plots/exampleSeqs/'];
end
%mkdir(out_folder);

if isfield(P(pt),'sz') == 0
    continue
end

if isfield(P(pt).sz(1),'seq_matrix') == 0 && isfield(P(pt),'seq_matrix') == 0
    continue
end



if example == 0
    dataName = P(pt).ieeg_name;
end
ptname = P(pt).name;


%% Concatenate all sequences
allSeq = [];
if isfield(P(pt),'seq_matrix') == 0
    if isempty(P(pt).sz(1).seq_matrix) == 1
        continue
    end
    for j = 1:length(P(pt).sz)
        allSeq = [allSeq,P(pt).sz(j).seq_matrix];
    end
else
    allSeq = P(pt).seq_matrix;
    if isempty(allSeq) == 1
        continue
    end
end


szTimes = P(pt).newSzTimes;

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
        if size(nonIctalSeq,2) < nseq 
            fprintf('Not enough sequences for %s, skipping\n',P(pt).name);
            continue
        end
        y = randsample(size(nonIctalSeq,2),nseq);
        seqs = nonIctalSeq(:,y);

    
    elseif ic == 1
        fprintf('Picking a random set of ictal sequences from all seizures\n');
        if size(ictalSeq,2) < nseq 
            fprintf('Not enough sequences for %s, skipping\n',P(pt).name);
            continue
        end
        y = randsample(size(ictalSeq,2),nseq);
        seqs = ictalSeq(:,y);
    end
else
    seqs = whichSeq;
    
end

%seqs = seqs_temp;

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
    if example == 0
        data = getiEEGData(dataName,0,0,pwfile);  
        fs = data.fs;
    

        %% Get the data for these times
        thresh.tmul = P(pt).thresh.tmul;
        thresh.absthresh = P(pt).thresh.absthresh;
        [gdf,extraoutput] = getSpikesSimple(P,pt,times,7,thresh,0);
        values = extraoutput.values;
    elseif example == 1
        gdf = P(pt).gdf;
        values = P(pt).eeg_data.values;
        fs = P(pt).fs;
        
        % Subsample the values to just be the times of interest
        example_times = P(pt).eeg_data.times;
        indices = max(1,round((times(1)-example_times(1))*fs)):...
            min(size(values,1),round((times(2)-example_times(1))*fs));
        values = values(indices,:);
    end
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
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m',...
    'b','r','g','c','m','b','r','g','c','m','b','r','g','c','m'};

if example == 0
for k = 1:length(seq_fig)
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

    
    for s = seq_fig{k}
        whichcol = mod(floor((s-1)/2)+1,5);
        whichcol(whichcol==0) = 5;
        whichrow = mod(s+1,2)+1;
        subplot('Position',[(whichcol-1)*1/columns (prows-whichrow)*1/prows 1/columns 1/prows])
        range = 0;
        pl = zeros(length(seq(s).whichCh),1);
        for i = 1:length(seq(s).whichCh)
            ch = seq(s).whichCh(i);
            amps = seq(s).values(:,ch) - mean(seq(s).values(:,ch)) - range;
            pl(i) = plot(seq(s).plottimes,amps,colors{i});
            hold on


            % find the times of the spikes with the desired channel
            spiketimes = seq(s).spikes(seq(s).spikes(:,1) == ch,2)-seq(s).time_col(1)+surroundtime;

            spikeamp = amps(round(fs*spiketimes));
            %spikeamp = ones(size(spiketimes,1),1)*max(seq(s).values(:,ch))-range;

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

        set(gca,'fontsize',10);
        xticks = 1:2:surroundtime*2-1;
        yl = ylim;
        yloc = yl(1)+(yl(2)-yl(1))*0.05;
        for t = 1:length(xticks)
           tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
        end

    end



    
if 0
    if save_plots == 1
        if isempty(whichSeq) == 0
            outputFile = [ptname,sprintf('_%d_%d_%d_%d-%d',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
        elseif isempty(whichSeq) == 1
            if ic == 1
                outputFile = [ptname,sprintf('_ic_%d_%d_%d_%d-%d',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
            elseif ic == 0
                outputFile = [ptname,sprintf('_inter_%d_%d_%d_%d-%d',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
            end
        end


        saveas(gcf,[out_folder,outputFile,'.png'])
        close(gcf)
    end
end

end
end


%% Plot again but make it pretty
for k = 1:length(seq_fig)
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

    for s = seq_fig{k}
        whichcol = mod(floor((s-1)/2)+1,5);
        whichcol(whichcol==0) = 5;
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

            spikeamp = amps(round(fs*spiketimes));

            %scatter(spiketimes,spikeamp,80,'k','filled');
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
    annotation('textbox',[0.45 0.905 0.1 0.1],'String',sprintf('Cluster %d',which_cluster),...
        'FitBoxToText','on','edgecolor','none','fontsize',20);




    if save_plots == 1
        if isempty(whichSeq) == 0
            outputFile = [ptname,sprintf('_%d_%d_%d_%d-%d_pretty',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
        elseif isempty(whichSeq) == 1
            if ic == 1
                outputFile = [ptname,sprintf('_ic_%d_%d_%d_%d-%d_pretty',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
            elseif ic == 0
                outputFile = [ptname,sprintf('_inter_%d_%d_%d_%d-%d_pretty',...
                    P(pt).thresh.whichDetector,...
                    P(pt).thresh.tmul,...
                    P(pt).thresh.absthresh,...
                    seq_fig{k}(1),seq_fig{k}(end))];
            end
        end    

        saveas(gcf,[out_folder,outputFile,'.png'])
        close(gcf)
    end

end



end

end