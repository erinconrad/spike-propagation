function bothSpikesAndSeqs(pt,whichPt,window)

% output file name
[~,~,~,resultsFolder,~] = fileLocations;



[chunk_seqs,times_plot_seq,MI,rl,dot,chunk_seqs_chs,vec,early] = seqFreqOverTime(pt,whichPt,window);
[chunk_spikes,times_plot_spike,MI_sp,chunk_spike_chs] = spikeFreqOverTime(pt,whichPt,window);

p_spike = zeros(length(chunk_spikes),1);
p_seq = zeros(length(chunk_seqs),1);
p_spike_all = [];
p_seq_all = [];

%% Spike freq
figure
for i = 1:size(chunk_spikes,2)
    times = times_plot_spike{i};
   p_spike(i) = plot(times/3600,chunk_spikes{i},'k','LineWidth',2);
    hold on  
    p_spike_all = [p_spike_all,p_spike(i)];
end

for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    p_seq(i) = plot(times/3600,chunk_seqs{i}*60,'k--','LineWidth',2);
    hold on  
    p_seq_all = [p_seq_all,p_seq(i)];
end




yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end

xlabel('Hour');
ylabel('Spike and sequence frequency');
legend([p_spike_all(1),p_seq_all(1),sz],{'Spikes per second','Sequences per minute','Seizure times'})
title(sprintf('Spike and sequence frequency over time for %s',pt(whichPt).name));
set(gca,'FontSize',15)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'Position',[50 100 1200 400])
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','SpikeSeqFreq/'];
mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'SpikeAndSeqFreq.png']);



%% MI
figure

for i = 1:size(chunk_spikes,2)
    times = times_plot_spike{i};
    sf = plot(times/3600,MI_sp{i},'k','LineWidth',2);
    hold on  
end

for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    so = plot(times/3600,MI{i},'k--','LineWidth',2);
    hold on  
end



yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end


xlabel('Hour');
ylabel('Moran index of spike frequency and recruitment latency');
legend([sf,so,sz],{'Spike frequency','Recruitment Latency','Seizure times'});
title(sprintf('Spatial autocorrelation of spike sequences over time for %s',pt(whichPt).name));
set(gca,'FontSize',15)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'Position',[50 100 1200 400])
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','MI/'];
mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'MI.png']);


%% Dot product
figure
for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    do = plot(times/3600,dot{i},'k','LineWidth',2);
    hold on  
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end

xlabel('Hour');
ylabel('Dot product');
legend([do,sz],{'Dot product','Seizure times'});
title(sprintf('Dot product of spike sequence vector with overall average for %s',pt(whichPt).name));
set(gca,'FontSize',15)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'Position',[50 100 1200 400])
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','dotProduct/'];
mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'dotProduct.png']);

%% Vector length
figure
for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    vl = plot(times/3600,vecnorm(vec{i},2,2),'k','LineWidth',2);
    hold on  
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end

xlabel('Hour');
ylabel('Vector length (mm)');
legend([vl,sz],{'Vector length (mm)','Seizure times'});
title(sprintf('Spike propagation vector length for %s',pt(whichPt).name));
set(gca,'FontSize',15)

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

set(gcf,'Position',[50 100 1200 400])
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','vecLength/'];
mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'vecLength.png']);


%% Make a couple movies (only for first time chunk)

saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','movies/'];
mkdir(saveFolder)
for i = 1:length(times_plot_seq{1})
info.title{i} = sprintf('Changing vectors for %s time %1.1f',pt(whichPt).name,times_plot_seq{1}(i));
end
info.save = [saveFolder,pt(whichPt).name,'vectors.gif'];
movieChangeVectors(early{1},vec{1},pt(whichPt).electrodeData.locs(:,2:4),info)

%{

% For each seizure, find which times are closest to the seizure
closeToSz = zeros(length(times_plot_spike{1}),1);
for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    szDist = abs(meanSzTimes - times_plot_spike{1});
    [~,I] = min(szDist);
    closeToSz(I) = 1;
end

for i = 1:length(times_plot_spike{1})
    if closeToSz(i) == 0
        title1{i} = sprintf('Spatial autocorrelation of\nspike frequency for time chunk %d',i);
    else
        title1{i} = ...
            sprintf('Spatial autocorrelation of\nspike frequency for time chunk %d\nSEIZURE',i);
    end
end
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','movies/'];
mkdir(saveFolder)
info.title = title1;
info.save = [saveFolder,pt(whichPt).name,'spike_MI.gif'];
brainMovieOfAnything(chunk_spike_chs{1},...
    pt(whichPt).electrodeData.locs(:,2:4),info)


for i = 1:length(times_plot_seq{1})
    if closeToSz(i) == 0
    title2{i} = sprintf('Spatial autocorrelation of\nrecruitment latency for time chunk %d',i);
    else
    title2{i} = sprintf('Spatial autocorrelation of\nrecruitment latency for time chunk %d\nSEIZURE',i);
    end
end
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','movies/'];
mkdir(saveFolder)
info.title = title2;
info.save = [saveFolder,pt(whichPt).name,'RL_MI.gif'];
brainMovieOfAnything(chunk_seqs_chs{1},...
    pt(whichPt).electrodeData.locs(:,2:4),info)

%}

end
