function aes_plot_aoi(allAllDist,allSADist,allFreqDist)

col1 = [0, 0.4470, 0.7410];
col2 = [0.8500 0.3250 0.0980];

n = length(allAllDist);
if size(allFreqDist,1) == 1
    allFreqDist = allFreqDist';
end
if size(allSADist,1) == 1
    allSADist = allSADist';
end

figure
xdata1 = ones(n,1).*rand(n,1)/2;
scatter(xdata1,allAllDist,150,col1,'filled')
hold on
plot([min(xdata1)-0.1 max(xdata1)+0.1],[median(allAllDist) median(allAllDist)],...
    'color',col1,'linewidth',3)
xdata2 = ones(n,1).*rand(n,1)/2+1*ones(n,1);
scatter(xdata2,allFreqDist,150,col2,'filled')
plot([min(xdata2)-0.1 max(xdata2)+0.1],[median(allFreqDist) median(allFreqDist)],...
    'color',col2,'linewidth',3)
plot([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2],...
    [max([allAllDist;allFreqDist])+10 max([allAllDist;allFreqDist])+10],...
    'k','linewidth',2)
text(((min(xdata1)+max(xdata1))/2+(min(xdata2)+max(xdata2))/2)/2-0.1,...
    max([allAllDist;allFreqDist])+13,'**','fontsize',50)
ylim([-10 max([allAllDist;allFreqDist])+20]);
xlim([-0.2 1.7])
xticks([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2]);
xticklabels({'All electrodes','Highest spike rate'})
ylabel('Distance from seizure onset zone (mm)')
set(gca,'fontsize',20)


figure
scatter(xdata1,allAllDist,150,col1,'filled')
hold on
plot([min(xdata1)-0.1 max(xdata1)+0.1],[median(allAllDist) median(allAllDist)],...
    'color',col1,'linewidth',3)
scatter(xdata2,allSADist,150,col2,'filled')
plot([min(xdata2)-0.1 max(xdata2)+0.1],[median(allSADist) median(allSADist)],...
    'color',col2,'linewidth',3)
plot([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2],...
    [max([allAllDist;allSADist])+10 max([allAllDist;allSADist])+10],...
    'k','linewidth',2)
text(((min(xdata1)+max(xdata1))/2+(min(xdata2)+max(xdata2))/2)/2-0.1,...
    max([allAllDist;allSADist])+13,'**','fontsize',50)
ylim([-10 max([allAllDist;allSADist])+20]);
xlim([-0.2 1.7])
xticks([(min(xdata1)+max(xdata1))/2,(min(xdata2)+max(xdata2))/2]);
xticklabels({'All electrodes','Largest area of influence'})
ylabel('Distance from seizure onset zone (mm)')
set(gca,'fontsize',20)

end