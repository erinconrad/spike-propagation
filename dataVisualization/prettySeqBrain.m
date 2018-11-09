function prettySeqBrain(pt,whichPt,whichSz,whichSeq)

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);

destFolder = [resultsFolder,'pretty_plots/Fig1/'];

if isempty(whichSz) ==1
    seq = pt(whichPt).seq_matrix(:,whichSeq);
else
    seq = pt(whichPt).sz(whichSz).seq_matrix(:,whichSeq);
end


% Get the non nan times for the sequence
spike_chs = find(isnan(seq) == 0);
spike_times = seq(find(isnan(seq) == 0));
[spike_times,I] = sort(spike_times);
spike_chs = spike_chs(I);

% Get channel locs
locs = pt(whichPt).electrodeData.locs(:,2:4);

%% Make plot
figure
set(gcf,'Position',[200 200 800 550]);
scatter3(locs(:,1),locs(:,2),locs(:,3),100,[0.6 0.6 0.6],'LineWidth',1.3);
hold on
scatter3(locs(spike_chs,1),locs(spike_chs,2),locs(spike_chs,3),100,'r','filled');

% Plot paths with arrows
for i = 1:length(spike_times)-1
    start = locs(spike_chs(i),:);
    finish = locs(spike_chs(i+1),:);
    dp = finish-start;
    quiver3([start(1)],[start(2)], [start(3)],...
        dp(1), dp(2), dp(3),'Color',[0 0 0],'LineWidth',2,'MaxHeadSize',1.5);
    if i == 1
        text(start(1),start(2)-12,start(3)-4,'Start','FontSize',25);
    end
    
end

finish = locs(spike_chs(end),:);
text(finish(1),finish(2)-5,finish(3)+5,'End','FontSize',25);
xticklabels([]);
yticklabels([]);
zticklabels([]);

view([0.7 0.2 0.2])

print(gcf,[destFolder,'brainSeq'],'-depsc');
eps2pdf([destFolder,'brainSeq.eps'])


%% Get vector
xyChan = pt(whichPt).electrodeData;
[vec,early_mean_all,late_mean_all] = getVectors2(seq,xyChan);




%% Make plot
figure
set(gcf,'Position',[200 200 800 550]);
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
%scatter3(locs(spike_chs,1),locs(spike_chs,2),locs(spike_chs,3),100,'r');
scatter3(locs(spike_chs(1),1),locs(spike_chs(1),2),locs(spike_chs(1),3),...
    100,'g','filled');

% Plot paths with arrows
%{
for i = 1:length(spike_times)-1
    start = locs(spike_chs(i),:);
    finish = locs(spike_chs(i+1),:);
    dp = finish-start;
    quiver3([start(1)],[start(2)], [start(3)],...
        dp(1), dp(2), dp(3),'Color',[0 0 0],'LineWidth',1,'MaxHeadSize',1); 
end
%}

quiver3(locs(spike_chs(1),1),locs(spike_chs(1),2),locs(spike_chs(1),3),...
    vec(1),vec(2),vec(3),'color',[0 0 0],'LineWidth',3,'MaxHeadSize',2); 


xticklabels([]);
yticklabels([]);
zticklabels([]);

view([0.7 0.2 0.2])

print(gcf,[destFolder,'brainVec'],'-depsc');
eps2pdf([destFolder,'brainVec.eps'])


end