function justPlotSeqs(pt,whichPts,window)

for whichPt = whichPts
% output file name
[~,~,~,resultsFolder,~] = fileLocations;



[chunk_seqs,times_plot_seq,MI,rl,angle,chunk_seqs_chs,vec,early] = seqFreqOverTime(pt,whichPt,window);
%[chunk_spikes,times_plot_spike,MI_sp,chunk_spike_chs] = spikeFreqOverTime(pt,whichPt,window);

%p_spike = zeros(length(chunk_spikes),1);
p_seq = zeros(length(chunk_seqs),1);
%p_spike_all = [];
p_seq_all = [];

%% Spike freq

figure


for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    p_seq(i) = plot(times/3600,chunk_seqs{i}*60,'k','LineWidth',2);
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
ylabel('Sequence frequency');
legend([p_seq_all(1),sz],{'Sequences per minute','Seizure times'})
title(sprintf('Sequence frequency over time for %s',pt(whichPt).name));
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
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','SeqFreq/'];
mkdir(saveFolder)
saveas(gcf,[saveFolder,pt(whichPt).name,'SeqFreq.png']);
close(gcf)


%% MI
figure

for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    so = plot(times/3600,MI{i},'k','LineWidth',2);
    hold on  
end



yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end


xlabel('Hour');
ylabel('Moran index of recruitment latency');
legend([so,sz],{'Recruitment Latency','Seizure times'});
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
close(gcf)

%% Dot product
figure
for i = 1:size(chunk_seqs,2)
    times = times_plot_seq{i};
    do = plot(times/3600,angle{i},'k','LineWidth',2);
    hold on  
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    sz = plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k-.','LineWidth',2);
end

xlabel('Hour');
ylabel('Angle (degrees)');
legend([do,sz],{'Angle','Seizure times'});
title(sprintf('Angle between spike sequence vector and overall average for %s',pt(whichPt).name));
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
close(gcf)

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
close(gcf)

%% Make a couple movies (only for first time chunk)

%{
saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','movies/'];
mkdir(saveFolder)
for i = 1:length(times_plot_seq{1})
info.title{i} = sprintf('Changing vectors for %s time %1.1f',pt(whichPt).name,times_plot_seq{1}(i));
end
info.save = [saveFolder,pt(whichPt).name,'vectors.gif'];
movieChangeVectors(early{1},vec{1},pt(whichPt).electrodeData.locs(:,2:4),info)



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
end
