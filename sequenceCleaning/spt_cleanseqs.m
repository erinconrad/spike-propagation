function [clueGraph,seq_track] = spt_cleanseqs(clueGraph)

% 'Noisy' sequences are detected using Kmeans clustering 
% and the weighted degree distribution. 

%Perform kmeans clustering
seq_track  = 1:size(clueGraph,2);                 %track removed sequences
in_deg     = sum(clueGraph,1);     
out_deg    = sum(clueGraph,2);     
degree     = in_deg + out_deg';    
K_data     = [in_deg',out_deg];                   %kmeans inputs
[idx,C]    = kmeans(K_data,3);                    %perform kmeans

%Which idx contains the noisy seqs?
grp1     = mean(degree(find(idx==1)));
grp2     = mean(degree(find(idx==2)));
grp3     = mean(degree(find(idx==3)));
grps     = [grp1,grp2,grp3];
[~,toss] = min(grps);
[~,highest] = max(grps);
middle = setdiff([1 2 3],[toss,highest]);

%Toss noisy sequences from clueGraph
toss_idx                     = find(idx==toss)';
clueGraph(:,toss_idx)        = [];
clueGraph(toss_idx,:)        = [];
seq_track(toss_idx)          = [];

%Plot an image of the sequence similarity matrix (clueGraph)
figure
imagesc(clueGraph)
colorbar
xlabel('Sequence number');
ylabel('Sequence number');
title('Sequence similarity matrix');
set(gca,'fontsize',15)
filen = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/results/SSM';
print(filen,'-depsc');
newPath = [filen,'.eps'];

%Plot a histogram of degree centrality and save it
figure
h = histogram(degree,round(length(idx)/3));
set(h,'facecolor','k');
hold on
lowborder = (max(degree(find(idx==toss)))+min(degree(find(idx==middle))))/2;
highborder = (max(degree(find(idx==middle)))+min(degree(find(idx==highest))))/2;
plot([lowborder lowborder],get(gca,'ylim'),'--','linewidth',3,'color','k');
plot([highborder highborder],get(gca,'ylim'),'--','linewidth',3,'color','k');
xlabel('Degree centrality');
ylabel('Count');
set(gca,'FontSize',15);

filen = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/results/histogram';
print(filen,'-depsc');
newPath = [filen,'.eps'];