function plotBrainPretty(saveFolder,szTimes,clusters,prop_c,colors,sparse_time,sparse_plot,sparse_cidx,plot_times,sum_times,locs,C,textLeg,pt,whichPt)

figure
set(gcf,'Position',[71 62 1084 717])
[ha, pos] = tight_subplot(3, 2, [.06 .01],[.1 .05],[.04 .01]);
delete(ha(2));
delete(ha(4));
delete(ha(6));
set(ha(1),'Position',[pos{1}(1) pos{1}(2) pos{1}(3)*1 pos{1}(4)]);
set(ha(3),'Position',[pos{3}(1) pos{3}(2) pos{3}(3)*2 pos{3}(4)]);
set(ha(5),'Position',[pos{5}(1) pos{5}(2) pos{5}(3)*2 pos{5}(4)]);

% Subplot 1: Plot the x, y, z over time
axes(ha(1));
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
hold on
plLeg = zeros(size(C,1),1);

for k = 1:size(C,1)
    plLeg(k) = scatter3(C(k,1),C(k,2),C(k,3),100,colors((k),:),'filled');
end
title(sprintf('Spike location centroids'),'fontsize',20)
set(gca,'FontSize',15);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
annotation('textbox',[0 0.77 0.2 0.2],'String','A','EdgeColor','none','fontsize',30);
legend(plLeg,textLeg,'Location','northeastoutside','fontsize',20);
view([55.3,2])


% Subplot 2: Plot the x, y, z over time
axes(ha(3));
ttext = {'x','y','z'};
toAdd = 0;






for i = 1:3
    % Plot each coordinate over time, offset from each other
    scatter(sparse_time/3600,sparse_plot(:,i)+...
        repmat(toAdd,size(sparse_plot,1),1),20,sparse_cidx)

    if i == 1
        minPoint = min(sparse_plot(:,i)+...
        repmat(toAdd,size(sparse_plot,1),1));
    end
    hold on
  %  text(0,toAdd+median(plot_thing(:,i)),...
  %      sprintf('%s',ttext{i}),'FontSize',30);
    tickloc(i) = toAdd+median(sparse_plot(:,i));
    if i == 3
        maxPoint = max(sparse_plot(:,i)+...
        repmat(toAdd,size(sparse_plot,1),1));
    end
    if i ~=3
        % Define the offset for each coordinate
        toAdd = toAdd + 300;
      %  toAdd = toAdd + 10+(max(sparse_plot(:,i)) - ...
      %      min(sparse_plot(:,i+1)));
    end
end
ylim([minPoint maxPoint + 30])
yticks(tickloc)
yticklabels({'X','Y','Z'});
set(gca,'fontsize',20);

% Plot the seizure times
%{
for j = 1:size(szTimes,1) 
    yl = ylim;
    plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
end
%}



set(gca,'xtick',[]);
xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
title(sprintf('X, Y, Z coordinates of all spikes'),'fontsize',20);







% Subplot 3: Plot the cluster distribution over time

axes(ha(5));
pl = zeros(length(clusters),1);
for i = 1:length(clusters)
    pl(i)= plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
        'color',colors((i),:),'LineWidth',2);
hold on
end

%{
for j = 1:size(szTimes,1) 
    yl = ylim;
    plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
end
%}
xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])

title(sprintf(['Proportion of spikes in given cluster '...
'(moving average)']));
xlabel('Time (hours)');

set(gca,'FontSize',20);
annotation('textbox',[0 0.5 0.2 0.2],'String','B','EdgeColor','none','fontsize',30);
annotation('textbox',[0 0.2 0.2 0.2],'String','C','EdgeColor','none','fontsize',30);

axes(ha(5));
% Plot the seizure times
for j = 1:size(szTimes,1) 
    yl = ylim;
    plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
        yl,'k','LineWidth',2);
end

axes(ha(3));
% Plot the seizure times
for j = 1:size(szTimes,1) 
    yl = ylim;
    plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
        yl,'k','LineWidth',2);
end

% Add annotations
annotation('textarrow',[0.2 0.275],[0.47 0.47],'String','Seizure',...
    'FontSize',20);

annotation('textarrow',[0.49-.075 0.495],[0.47 0.47],'String','Seizure',...
    'FontSize',20);

annotation('textarrow',[0.735-.075 0.74],[0.47 0.47],'String','Seizure',...
    'FontSize',20);



%pause

print(gcf,[saveFolder,'clustTimePretty_',sprintf('%s',pt(whichPt).name)],'-depsc');
eps2pdf([saveFolder,'clustTimePretty_',sprintf('%s',pt(whichPt).name),'.eps'])
close(gcf)

end