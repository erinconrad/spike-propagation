function CPlots(pt,cluster,whichPts)

%{ 

CPlots
This is my cleaned up file for getting plots on the cluster data 

%}    

doPretty = 1;
makeSparse = 0;

% Save file location
[~,~,~,resultsFolder,~] = fileLocations;
destFolder = [resultsFolder,'clustering/plots/'];
mkdir(destFolder);


if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0 
            whichPts = [whichPts,i];
        end
    end
elseif whichPts == 100
    whichPts = [1,4 6 8 9 15 17 18 19 20 22 24 25 27 30 31];
end

allCounts = [];
allPat = [];
allChunk = [];


for whichPt = whichPts
    
    fprintf('Doing %s\n',pt(whichPt).name);
    
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder);
    
    % Get patient parameters
    locs = pt(whichPt).electrodeData.locs(:,2:4);
    szTimes = pt(whichPt).newSzTimes;
    soz = pt(whichPt).newSOZChs;
    
    % Reorder seizure times if out of order
    oldSzTimes = szTimes;
    szTimes = sort(szTimes,1);
    if isequal(oldSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s seizure times out of order\n',pt(whichPt).name);
    end
    
    % Combine nearly equal seizure times
    newIdx = 2;
    newSzTimes = [];
    newSzTimes(1,:) = szTimes(1,:);
    for j = 2:size(szTimes,1)
        if abs(szTimes(j,1)-szTimes(j-1,1)) < 10 && ...
                abs(szTimes(j,2)-szTimes(j-1,2))
           newIdx = newIdx - 1; 
           newSzTimes(newIdx,1) = min(szTimes(j,1),szTimes(j-1,1));
           newSzTimes(newIdx,2) = max(szTimes(j,2),szTimes(j-1,2));  
        else
           newSzTimes(newIdx,:) = szTimes(j,:);
        end
        newIdx = newIdx + 1;
    end
    
    if isequal(newSzTimes,szTimes) == 0
        fprintf('WARNING!!! %s had duplicate seizure times\n',pt(whichPt).name);
    end
    
    % Pull cluster info
    %{
    all_times_all = pt(whichPt).cluster.all_times_all; % all spike times
    all_spikes = pt(whichPt).cluster.all_spikes; % all spike channels
    all_locs = pt(whichPt).cluster.all_locs;
    k = pt(whichPt).cluster.k; % the number of clusters
    idx = pt(whichPt).cluster.idx; % the cluster index for every spike
    C = pt(whichPt).cluster.C; % the centroids of the clusters
    bad_cluster = pt(whichPt).cluster.bad_cluster; % which clusters are bad
    %}
    all_times_all = cluster(whichPt).all_times_all; % all spike times
    all_spikes = cluster(whichPt).all_spikes; % all spike channels
    all_locs = cluster(whichPt).all_locs;
    k = cluster(whichPt).k; % the number of clusters
    idx = cluster(whichPt).idx; % the cluster index for every spike
    C = cluster(whichPt).C; % the centroids of the clusters
    bad_cluster = cluster(whichPt).bad_cluster; % which clusters are bad
    
    % Confirm that I do not have any ictal spikes
    t = find(any(all_times_all >= szTimes(:,1)' & all_times_all <= szTimes(:,2)',2));
    if isempty(t) == 0
        fprintf('WARNING: Remaining ictal spikes for %s!\n',pt(whichPt).name);
        all_times_all(t) = [];
        all_spikes(t) = [];
        all_locs(t,:) = [];
        idx(t) = [];
    end
    
    % Skip patient if all clusters are bad
    if isequal(1:k,bad_cluster) == 1
        fprintf('All clusters for %s are bad, skipping.\n',pt(whichPt).name);
        continue
    end

    % Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    C(bad_cluster,:) = [];
    
    %% Do plots
    %colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
    %colors = [0 0 1;1 0 0;0.7 0.1840 0.6];
    colors = [0 0 1;1 0 0;0.9290 0.540 0.1250];
    
    c_idx = zeros(size(idx,1),3);
    
    
    
    if doPretty == 1
        
        
        plot_thing = all_locs;
        plot_times = all_times_all;
        % Get the times for spikes in each cluster
        for i = 1:length(clusters)
            clust{i} = plot_times(idx == clusters(i));
        end
        window = 3600;
        
        
        [sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
        totalSum = zeros(1,size(sum_times,2));
        for i = 1:length(clusters)
            totalSum = totalSum + sum_c(i,:);
        end
        prop_c = sum_c./totalSum;
        %}
        
        for i = 1:length(clusters)  
           c_idx(idx==clusters(i),:) = ...
               repmat(colors(i,:),sum(idx==clusters(i)),1);
        end
        
        possibleText = {'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6'};
        textLeg = possibleText(1:size(C,1));

        
        %% Figure 
        
        sparse_plot = plot_thing;
        sparse_time = plot_times-min(plot_times);
        sparse_cidx = c_idx;
        %plotBrainPretty(saveFolder,szTimes,clusters,prop_c,colors,sparse_time,sparse_plot,sparse_cidx,plot_times,sum_times,locs,C,textLeg,pt,whichPt);
        %plotHourBins(szTimes,clusters,prop_c,colors,sparse_time,plot_times,sum_times);
        
        figure
        set(gcf,'position',[1 409 1440 320]);
        for i = 1:length(clusters)
            pl(i)= plot((sum_times-min(plot_times))/3600,prop_c(i,:),...
                'color',colors((i),:),'LineWidth',2);
        hold on
        end

        
        for j = 1:size(szTimes,1) 
            yl = ylim;
          szp = plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
                yl,'k--','LineWidth',3);
        end
        
        xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
        %{
        title(sprintf(['Proportion of spikes in given cluster, '...
        'moving average']));
        %}
        xlabel('Time (hours)');
        %{
        title(sprintf(['Proportion of spikes in given cluster '...
        '(moving average)']))
        %}
        set(gca,'FontSize',25);
        
        if whichPt == 9
            legend([pl,szp],[textLeg,'Seizures'],'Position',[0.65 0.35 0.15 0.1])
            
        elseif whichPt == 31
            legend([pl,szp],[textLeg,'Seizures'],'Position',[0.83 0.36 0.1 0.1])
        elseif whichPt == 17
            legend([pl,szp],[textLeg,'Seizures'],'Position',[0.4 0.5 0.1 0.1])
        end
        
        print(gcf,[saveFolder,'new_pretty_',sprintf('%s',pt(whichPt).name)],'-depsc');
        
        
        figure
        set(gcf,'Position',[71 241 1237 538])
        [ha, ~] = tight_subplot(2, 1, [.09 .03],[.14 .07],[.05 .01]);
        
        % Subplot 1: Plot the x, y, z over time
        axes(ha(1));
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
                if whichPt == 9
                    toAdd = toAdd + 100;
                elseif whichPt == 31
                    toAdd = toAdd + 270;
                end
                if i == 2 && whichPt == 9
                    toAdd = toAdd + 300;
                end
                if whichPt ~= 9 && whichPt ~= 31
                    toAdd = toAdd + 10+(max(sparse_plot(:,i)) - ...
                   min(sparse_plot(:,i+1)));
                end
            end
        end
        ylim([minPoint maxPoint + 30])
     %   yticks(tickloc)
        yticklabels({'X','Y','Z'});
        title(sprintf('X, Y, Z coordinates of all spikes for %s',...
            pt(whichPt).name));
        set(gca,'fontsize',25);

        % Plot the seizure times
        %{
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
        end
        %}
        
        

        set(gca,'xtick',[]);
        xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
        %{
        title(sprintf('X, Y, Z coordinates of all spikes for %s',...
            pt(whichPt).name),'fontsize',20);
        %}
        
        
        

       % annotation('textbox',[0 0.8 0.2 0.2],'String','D','EdgeColor','none','fontsize',30);

        
        % Subplot 2: Plot the cluster distribution over time
        
        axes(ha(2));
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
        %{
        title(sprintf(['Proportion of spikes in given cluster, '...
        'moving average']));
        %}
        xlabel('Time (hours)');
        title(sprintf(['Proportion of spikes in given cluster '...
        '(moving average)']))
    
        
       
        set(gca,'FontSize',25);
       % annotation('textbox',[0 0.347 0.2 0.2],'String','E','EdgeColor','none','fontsize',30);
        
        axes(ha(1));
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
                yl,'k','LineWidth',2);
        end
        
        axes(ha(2));
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,...
                yl,'k','LineWidth',2);
        end
        
        if whichPt == 9
            legend(pl,textLeg,'Position',[0.65 0.35 0.15 0.1])
            
        elseif whichPt == 31
            legend(pl,textLeg,'Position',[0.83 0.36 0.1 0.1])
        end
        
        % Add annotations
        if whichPt == 31
            annotation('textarrow',[0.22 0.285],[0.69 0.69],'String','Seizure',...
                'FontSize',25);

            annotation('textarrow',[0.5-.075 0.5],[0.69 0.69],'String','Seizure',...
                'FontSize',25);

            annotation('textarrow',[0.74-.075 0.744],[0.69 0.69],'String','Seizure',...
                'FontSize',25);
        elseif whichPt == 9
            annotation('textarrow',[0.08+0.075 0.08],[0.69 0.69],'String','Seizure',...
                'FontSize',25);

            annotation('textarrow',[0.5-.075 0.5],[0.69 0.69],'String','2 Seizures',...
                'FontSize',25);

            annotation('textarrow',[0.955-.075 0.955],[0.69 0.69],'String','Seizure',...
                'FontSize',25);
        elseif whichPt == 20
        end
        
        %{
        annotation('textbox',[0.32 0.95 0.5 0.05],'string',...
            sprintf('X, Y, Z coordinates of all spikes for %s',...
            pt(whichPt).name),'fontsize',23,'edgecolor','none');
        
        annotation('textbox',[0.29 0.49 0.5 0.05],'string',...
            sprintf(['Proportion of spikes in given cluster, '...
        'moving average']),'fontsize',23,'edgecolor','none');
            %|
        %}
        
        pause

        print(gcf,[saveFolder,'clustTimePretty_',sprintf('%s',pt(whichPt).name)],'-depsc');
        eps2pdf([saveFolder,'clustTimePretty_',sprintf('%s',pt(whichPt).name),'.eps'])
        close(gcf)
        

        
        
        
    else
        % Assign the sequence a color based on its cluster index
        
        for i = 1:length(idx)
           c_idx(i,:) = colors(idx(i),:); 
        end

        figure
        set(gcf,'Position',[50 100 1200 1200])

        %% Subplot 1: Plot the x, y, z over time
        subplot(3,1,1)
        ttext = {'x','y','z'};
        toAdd = 0;
        plot_thing = all_locs;
        plot_times = all_times_all;

        for i = 1:3
            % Plot each coordinate over time, offset from each other
            scatter(plot_times/3600,plot_thing(:,i)+...
                repmat(toAdd,size(plot_thing,1),1),20,c_idx)
            hold on
            text(plot_times(1)/3600-0.3,toAdd+median(plot_thing(:,i)),...
                sprintf('%s',ttext{i}),'FontSize',30);
            if i ~=3
                % Define the offset for each coordinate
                toAdd = toAdd + 10+(max(plot_thing(:,i)) - ...
                    min(plot_thing(:,i+1)));
            end
        end

        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
        end

        set(gca,'ytick',[]);
        xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
        title(sprintf('X, y, z coordinates of all spikes for %s',...
            pt(whichPt).name));
        set(gca,'FontSize',15);

        %% Subplot 2: Plot the cluster distribution over time

        % Get the times for spikes in each cluster
        for i = 1:length(clusters)
            clust{i} = plot_times(idx == clusters(i));
        end


        window = 1800;

        % NEED TO CHECK THIS SCRIPT - ALSO IT's SUUUPER SLOW
        
        [sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
        totalSum = zeros(1,size(sum_times,2));
        for i = 1:length(clusters)
            totalSum = totalSum + sum_c(i,:);
        end
        prop_c = sum_c./totalSum;
        %}


        % alternate method - try non-moving window

        % enough bins so that it's essentially 30 minute windows
        %{
        nbins = round((max(plot_times)-min(plot_times))/(window));
        [Y,E] = discretize(plot_times,nbins);
        new_times = E(2:end);
        new_counts = zeros(nbins,length(clusters));
        for bb = 1:nbins
            for k = 1:length(clusters)
                new_counts(bb,k) = sum(Y==bb & idx==clusters(k));
            end
        end
        new_prop = new_counts./sum(new_counts,2);
        %sum_times = new_times;
        %prop_c = new_prop';

        subplot(3,1,2)
        pl = zeros(length(clusters),1);
        for i = 1:length(clusters)
            pl(i)= plot(sum_times/3600,prop_c(i,:),...
                'color',colors(clusters(i),:),'LineWidth',2);
        hold on
        end

        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1) szTimes(j,1)]/3600,yl,'k','LineWidth',2);
        end
        xlim([plot_times(1)/3600-1 plot_times(end)/3600+1])
        title(sprintf(['Proportion of sequences in given cluster, '...
        '%d minute bins, %s'],window/60,pt(whichPt).name));
        set(gca,'FontSize',15);
        %}


        %% Subplot 3: Plot locations of centroids
        %{
        subplot(3,1,3)
        scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k');
        hold on
        for k = 1:size(C,1)
            scatter3(C(k,1),C(k,2),C(k,3),60,colors(clusters(k),:),'filled');
        end
        title(sprintf('Spike location centroids for each cluster for %s',...
            pt(whichPt).name))
        set(gca,'FontSize',15);
        set(gca,'xticklabel',[])
        set(gca,'yticklabel',[])
        set(gca,'zticklabel',[])


        %pause

        print(gcf,[saveFolder,'clustTime_',sprintf('%s',pt(whichPt).name)],'-depsc');
        eps2pdf([saveFolder,'clustTime_',sprintf('%s',pt(whichPt).name),'.eps'])
        close(gcf)
        %}
    end
end

end