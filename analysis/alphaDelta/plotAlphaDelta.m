function plotAlphaDelta(pt,cluster,power,whichPts)

%{

Analysis: does the alpha delta ratio correlate with cluster distribution?

Spearman Rank correlation - alpha/delta ratio in a 2000 s bin
against proportion of spikes in most popular cluster

%}



doPretty = 1;
skipDone = 1;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
destFolder = [resultsFolder,'alphaDelta/plots/'];
mkdir(destFolder);

if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
elseif whichPts == 100
    whichPts = [4,6,8,9,12,15,17,18,19,20,22,24,25,27,30,31];
end

for whichPt = whichPts
    
    %look = mean(power(whichPt).ad_rat,2);
    
    
    saveFolder = [destFolder,pt(whichPt).name,'/'];
    mkdir(saveFolder)
    fprintf('Doing %s\n',pt(whichPt).name);
    szTimes = pt(whichPt).newSzTimes;
    locs = pt(whichPt).electrodeData.locs(:,2:4);

    
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
    

    % Remove bad clusters
    bad_idx = find(ismember(idx,bad_cluster));
    all_times_all(bad_idx) = [];
    all_spikes(bad_idx) = [];
    all_locs(bad_idx,:) = [];
    idx(bad_idx) = [];
    clusters = 1:k; clusters(bad_cluster) = [];
    
    % Get most popular cluster
    popular = mode(idx);
    
    bin_times = pt(whichPt).runTimes;
    prop_pop = zeros(size(bin_times,1),1);
    
    % Run through bin times and get proportion of spikes in most popular
    % cluster for that bin.
    for i = 1:size(bin_times,1)
        
        % Cluster identities of spikes in between those times
        whichClust = idx(all_times_all > bin_times(i,1) & ...
            all_times_all < bin_times(i,2));
        
        prop_pop(i) = sum(whichClust == popular)/length(whichClust);
        
    end
    
    
    %[badChNums,badChNamesOut] = getBadChs(pt,whichPt);
    %mean_ad = mean(power(whichPt).ad_rat,1)';
    all_ad = power(whichPt).alpha./power(whichPt).delta;
    mean_ad = mean(all_ad,1);
    %{
    mean_ad = mean(power(whichPt).ad_rat(...
        ~ismember(1:length(pt(whichPt).channels),badChNums),:),1)';
    %}
    
    if 1 == 0
        figure
        subplot(2,1,1)
        plot(power(whichPt).times/3600,prop_pop);
        hold on
        plot(power(whichPt).times/3600,mean_ad);


        subplot(2,1,2)
        scatter(mean_ad,prop_pop)
    end
    
    mean_ad(isnan(prop_pop)) = [];
    prop_pop(isnan(prop_pop)) = [];
    
    
    [rho,pval] = corr(mean_ad',prop_pop,'Type','Spearman');
    fprintf(['For %s, the correlation between proportion in most popular'...
        ' cluster and alpha delta ratio is:\n %1.1f (p = %1.1e)\n'],...
        pt(whichPt).name,rho,pval);
    
    
    colors = [0 0 1;1 0 0;0 1 0; 0.5 0.5 1; 1 0.5 0.5; 0.5 1 0.5; 0.4 0.7 0.4];
    c_idx = zeros(size(idx,1),3);
    
    
    
    if doPretty == 1
        
        %% Prep for plot
        
        plot_thing = all_locs;
        plot_times = all_times_all;
        
        window = 3600;
        
        for i = 1:length(clusters)
            clust{i} = plot_times(idx == clusters(i));
        end
        
        [sum_c,sum_times] = movingSumCounts(clust,plot_times,window);
        totalSum = zeros(1,size(sum_times,2));
        for i = 1:length(clusters)
            totalSum = totalSum + sum_c(i,:);
        end
        prop_c = sum_c./totalSum;
        
        %% Plot
        
        figure
        set(gcf,'Position',[50 100 1200 500])
        [ha, ~] = tight_subplot(3, 1, [.06 .01],[.12 .08],[.05 .01]);
        
        % Subplot 1: Plot the x, y, z over time
        axes(ha(1));
        ttext = {'x','y','z'};
        toAdd = 0;
        
        
        
        for i = 1:length(clusters)  
           c_idx(idx==clusters(i),:) = ...
               repmat(colors(i,:),sum(idx==clusters(i)),1);
        end
        
        
        
        sparse_plot = plot_thing;
        sparse_time = plot_times-min(plot_times);
        sparse_cidx = c_idx;
       
        
        
        for i = 1:3
            % Plot each coordinate over time, offset from each other
            scatter(sparse_time/3600,sparse_plot(:,i)+...
                repmat(toAdd,size(sparse_plot,1),1),20,sparse_cidx)
            hold on
          %  text(0,toAdd+median(plot_thing(:,i)),...
          %      sprintf('%s',ttext{i}),'FontSize',30);
            tickloc(i) = toAdd+median(sparse_plot(:,i));
            if i ~=3
                % Define the offset for each coordinate
                toAdd = toAdd + 10+(max(sparse_plot(:,i)) - ...
                    min(sparse_plot(:,i+1)));
            end
        end
        yticks(tickloc)
        yticklabels({'X','Y','Z'});

        % Plot the seizure times
        %{
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,yl,'k','LineWidth',2);
        end
        %}
        
        ylim([min(ylim) max(ylim)+70])

        set(gca,'xtick',[]);
        xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
        title(sprintf('X, Y, Z coordinates of all spikes for %s',...
            pt(whichPt).name));
        
        
        
        set(gca,'FontSize',20);
        annotation('textbox',[0 0.79 0.2 0.2],'String','A','EdgeColor','none','fontsize',30);
        
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
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,yl,'k','LineWidth',2);
        end
        %}
        xlim([sparse_time(1)/3600-1 sparse_time(end)/3600+1])
        title(sprintf(['Proportion of spikes in given cluster, '...
        'moving average']));
        set(gca,'xtick',[]);
        %xlabel('Time (hours)');
        
        set(gca,'FontSize',20);
        annotation('textbox',[0 0.5 0.2 0.2],'String','B','EdgeColor','none','fontsize',30);
        
        
        axes(ha(3));
        plot((power(whichPt).times-min(plot_times))/3600,mean(all_ad,1),...
            'k','LineWidth',2);
        xlim([(plot_times(1)-min(plot_times))/3600-1 (plot_times(end)-min(plot_times))/3600+1])
        hold on
        
        title(sprintf('Alpha-delta power ratio averaged across all electrodes'))
        xlabel('Time (hours)');
        set(gca,'FontSize',20);
        annotation('textbox',[0 0.18 0.2 0.2],'String','C','EdgeColor','none','fontsize',30);

        
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,yl,'k','LineWidth',2);
        end
        
        axes(ha(1));
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,yl,'k','LineWidth',2);
        end
        
        axes(ha(2));
        % Plot the seizure times
        for j = 1:size(szTimes,1) 
            yl = ylim;
            plot([szTimes(j,1)-min(plot_times) szTimes(j,1)-min(plot_times)]/3600,yl,'k','LineWidth',2);
        end
        
        
        annotation('textarrow',[0.24 0.19],[0.92 0.92],'String','Seizures',...
            'FontSize',20);
        annotation('textarrow',[0.84 0.89],[0.92 0.92],'String','Seizure',...
            'FontSize',20);
        
        
        %pause

        print(gcf,[saveFolder,'ad_',sprintf('%s',pt(whichPt).name)],'-depsc');
        eps2pdf([saveFolder,'ad_',sprintf('%s',pt(whichPt).name),'.eps'])
        close(gcf)
        
    end

        
        
    
end

%figure
%plot(power(whichPt).times/3600,mean(power(whichPt).ad_rat,1));


end