function plotHourBins(szTimes,clusters,prop_c,colors,sparse_time,plot_times,sum_times)
   
    figure
    set(gcf,'Position',[440 457 1000 341]);
    for i = 1:length(clusters)
        plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
            'color',colors((i),:),'LineWidth',3);
        hold on
    end

    %{
    for j = 1:size(szTimes,1) 
        yl = ylim;
        plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
    end
    %}
    xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
    %{
    title(sprintf(['Proportion of spikes in given cluster, '...
    'moving average']));
    %}
    xlabel('Time (hours)');
    ylabel({'Proportion of spikes in cluster'});
    title('Proportion of spikes in each cluster (moving average)')

    set(gca,'FontSize',20);

    n_hours = floor(sparse_time(end)/3600+1);
    for i = 1:n_hours
        plot([i i],get(gca,'ylim'),'k--','Linewidth',1);
    end
    print('hour_bin_78','-depsc')


end