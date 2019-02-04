function prettyClusterSeqs(pt,cluster)


whichPt =30;
n_seqs = 4;
surroundtime = 2;


%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);


seq_matrix = pt(whichPt).seq_matrix;
seq_index = cluster(whichPt).seq_index;
idx = cluster(whichPt).idx;
k = cluster(whichPt).k;

% Remove ties
keep = ones(size(seq_matrix,2),1);
for s = 1:size(seq_matrix,2)
   curr_seq = seq_matrix(:,s);
   nonans = curr_seq(~isnan(curr_seq));
   norepeats = unique(nonans);
   if length(norepeats) < 0.5*length(nonans)
       keep(s) = 0;
   end
end
seq_matrix(:,keep==0) = [];


rep_seq = [];
for i = 1:k
    % the indices of the spikes in this cluster
    spike_clust = find(idx==i); 

    % Take 50 of these indices randomly
    whichSpikes = spike_clust(randperm(length(spike_clust),...
        min(n_seqs,length(spike_clust))));

    % Get the sequences these spikes belong to
    whichSeqs = seq_index(whichSpikes);
    rep_seq{i} = seq_matrix(:,whichSeqs);
    
end

%% Define the structure that I will put plotting data into
nseq = struct;

for i = 1:k
    for j = 1:n_seqs
        s = rep_seq{i}(:,j);
        
        % Get the non nan times for the sequence
        spike_chs = find(isnan(s) == 0);
        spike_times = s(find(isnan(s) == 0));

        % sort them from earliest to latest times in the sequence
        [spike_times,I] = sort(spike_times);
        spike_chs = spike_chs(I);

        times = [spike_times(1)-surroundtime,spike_times(1)+surroundtime];
        
        
        data = getiEEGData(pt(whichPt).ieeg_name,0,0,pwfile);  
        fs = data.fs;

        %% Get the data for these times
        thresh.tmul = pt(whichPt).thresh.tmul;
        thresh.absthresh = pt(whichPt).thresh.absthresh;
        [~,extraoutput] = getSpikesSimple(pt,whichPt,times,7,thresh,0);
        values = extraoutput.values;
        plottimes =  [1:size(values,1)]/fs;
        
        %% Get the spike times
        spikes = [spike_chs,spike_times];
        
        %% Load into the struct
        nseq(i).seq(j).seq = s;
        nseq(i).seq(j).spikes = spikes;
        nseq(i).seq(j).plottimes = plottimes;
        nseq(i).seq(j).values = values;
        nseq(i).seq(j).whichCh = spike_chs;
        nseq(i).seq(j).time_col = spike_times;
        
    end
end

figure
set(gcf,'Position',[1 248 1433 550]);
[ha, pos] = tight_subplot(k, n_seqs, 0.01, [0.05 0.05], [0.05 0.01]);
col = [0 0 1;1 0 0;0 1 0];
for j = 1:k
    seq = nseq(j).seq;
    for s = 1:n_seqs
        sp = (j-1)*n_seqs + s;
        axes(ha(sp))
        
        range = 0;
        pl = zeros(length(seq(s).whichCh),1);
        for i = 1:length(seq(s).whichCh)
            ch = seq(s).whichCh(i);
            amps = seq(s).values(:,ch) - mean(seq(s).values(:,ch)) - range;
            pl(i) = plot(seq(s).plottimes,amps,'color','k','linewidth',2);
            hold on

            % find the times of the spikes with the desired channel
            spiketimes = seq(s).spikes(seq(s).spikes(:,1) == ch,2)-seq(s).time_col(1)+surroundtime;

            spikeamp = amps(round(fs*spiketimes));
            %spikeamp = ones(size(spiketimes,1),1)*max(seq(s).values(:,ch))-range;

            scatter(spiketimes,spikeamp,80,'k','filled');
            range = range + max(seq(s).values(:,ch)) - min(seq(s).values(:,ch));   
        end
        if j == k
            xlabel('Time (s)');
        else
            set(gca,'XTickLabel',[]);
        end
        
        if s == 1
            set(gca,'ytick',mean(ylim));
            set(gca,'yticklabel',sprintf('Cluster %d',j));
        else
            set(gca,'YTickLabel',[]);
        end
        
        set(gca,'FontSize',15);
        
    end 
end
annotation('textbox',[0.43 0.87 0.1 0.1],'String','Example sequences',...
    'FitBoxToText','on','FontSize',20,'LineStyle','none');

print(gcf,'ex_seqs','-depsc')

end