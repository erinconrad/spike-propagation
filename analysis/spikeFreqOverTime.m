function spikeFreqOverTime(pt,whichPt,window)

% Get spike frequency overall and by channel

%% Remove EKG artifact and depth electrodes
rmEKG = 1;
rmDepth = 1;

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
gdfFolder = [resultsFolder,'gdf/'];


nchs = length(pt(whichPt).channels);

gdf_all{1} = [];


%% fill up the info from the first seizure
for k = 1:length(pt(whichPt).sz(1).chunkFiles)
        
    
    % Load gdf file
    load([gdfFolder,pt(whichPt).name,'/',pt(whichPt).sz(1).chunkFiles{k}]);

    if isempty(gdf) == 1
        continue
    end

    % Load gdf ekg file
    load([gdfFolder,pt(whichPt).name,'/',pt(whichPt).sz(1).EKGchunkFiles{k}]);

    % remove EKG artifact and depth electrodes
    if rmEKG == 1
        %gdf = removeEKGArtifact(gdf,gdf_ekg,prox);
    end

    if rmDepth == 1
        %gdf = removeChs(gdf,electrodeData,rmType);
    end
    
    gdf_all{1} = [gdf_all{1};gdf];

end

if length(pt(whichPt).sz) > 1
for j = 2:length(pt(whichPt).sz)
    
    gdf_curr_sz = [];
    
    %% Define time windows
    totalTime = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j).runTimes(1,1);
   
    szTime = pt(whichPt).sz(j).onset;
    szTimeLast = pt(whichPt).sz(j-1).onset;
    
    
    if isfield(pt(whichPt).sz,'runTimes') == 0
        continue
    end
    
    for k = 1:length(pt(whichPt).sz(j).chunkFiles)
        
        if exist([gdfFolder,pt(whichPt).name,'/',pt(whichPt).sz(j).chunkFiles{k}],'file') == 0
            continue
        end
        
        % Load gdf file
        load([gdfFolder,pt(whichPt).name,'/',pt(whichPt).sz(j).chunkFiles{k}]);

        if isempty(gdf) == 1
            continue
        end

        % Load gdf ekg file
        load([gdfFolder,pt(whichPt).name,'/',pt(whichPt).sz(j).EKGchunkFiles{k}]);

        % remove EKG artifact
        if rmEKG == 1
          %  gdf = removeEKGArtifact(gdf,gdf_ekg,prox);
        end
        
        if rmDepth == 1
          %  gdf = removeChs(gdf,electrodeData,rmType);
        end
        
        gdf_curr_sz = [gdf_curr_sz;gdf];
        
    end
   
    if szTime - szTimeLast > totalTime
        
        % If the seizure time is more than 24 hours after the last seizure,
        % put this in a new chunk for plotting purposes
        gdf_all{end+1} = gdf_curr_sz;
    else
        % if it is less than 12 hours after the last seizure, then need to
        % add any spikes that occured after the last seizure run time to
        % this new chunk
        keepAfter = pt(whichPt).sz(j-1).runTimes(end,2);
        gdf_to_keep = gdf_curr_sz(gdf_curr_sz(:,2) > keepAfter,:);
        gdf_all{end} = [gdf_all{end};gdf_to_keep];
    end
      
end
end

chunk_spikes = cell(size(gdf_all));
chunk_spike_chs = cell(size(gdf_all));
times_plot = cell(size(gdf_all));

%% Now divide spikes into windows
for i = 1:size(gdf_all,2)
    gdf = gdf_all{i};
    totalTime = gdf(end,2)-gdf(1,2);
    nchunks = ceil(totalTime/window);
    
    chunk_spikes{i} = zeros(nchunks,1);
    chunk_spike_chs{i} = zeros(nchunks,nchs); 
    times_plot{i} = zeros(nchunks,1);
    
    for tt = 1:nchunks
        times = [(tt-1)*window + gdf(1,2),tt*window + gdf(1,2)];
        times_plot{i}(tt) = (times(1)+times(2))/2;
        
        % get the appropriate spikes in this time
         correct_spikes =  gdf(gdf(:,2) >= times(1) & gdf(:,2) <= times(2),:);
         chunk_spikes{i}(tt) = size(correct_spikes,1);
         
         % get the number of spikes per channel
         for k = 1:size(correct_spikes,1)
            ch = correct_spikes(k,1);
            chunk_spike_chs{i}(tt,ch) = chunk_spike_chs{i}(tt,ch) + 1;
             
         end
        
    end
    chunk_spikes{i} = chunk_spikes{i}/window;
    chunk_spike_chs{i} = chunk_spike_chs{i}/window;
    
end

%% Plot
figure
for i = 1:size(chunk_spikes,2)
    times = times_plot{i};
    plot(times/3600,chunk_spikes{i},'k');
    hold on  
end

yl = ylim;

for j = 1:length(pt(whichPt).sz)
    szTimes = [pt(whichPt).sz(j).onset,pt(whichPt).sz(j).offset];
    meanSzTimes = (szTimes(1) + szTimes(2))/2;
    plot([meanSzTimes meanSzTimes]/3600,[yl(1) yl(2)],'k--');
end

xlabel('Hour');
ylabel('Spikes per hour');
title('Spike frequency over time');
set(gca,'FontSize',15)

end