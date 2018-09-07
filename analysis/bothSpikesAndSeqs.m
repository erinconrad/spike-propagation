function bothSpikesAndSeqs(pt,whichPt,window)

% output file name
[~,~,~,resultsFolder,~] = fileLocations;

saveFolder = [resultsFolder,'plots/',pt(whichPt).name,'/','SpikeSeqFreq/'];
mkdir(saveFolder)
saveFile = [saveFolder,pt(whichPt).name,'SpikeAndSeqFreq.png'];


[chunk_seqs,times_plot_seq] = seqFreqOverTime(pt,whichPt,window);
[chunk_spikes,times_plot_spike] = spikeFreqOverTime(pt,whichPt,window);

p_spike = zeros(length(chunk_spikes),1);
p_seq = zeros(length(chunk_seqs),1);
p_spike_all = [];
p_seq_all = [];

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
legend([p_spike_all,p_seq_all,sz],{'Spikes per second','Sequences per minute','Seizure times'})
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

saveas(gcf,saveFile);


end