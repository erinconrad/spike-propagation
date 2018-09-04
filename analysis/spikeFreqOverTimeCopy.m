function spikeFreqOverTime(pt,whichPt,window)

% Get spike frequency overall and by channel

%% Remove EKG artifact and depth electrodes
rmEKG = 1;
rmDepth = 1;


nchs = length(pt(whichPt).channels);

spike_freq = cell(length(pt(whichPt).sz),1);
spike_freq_ch = cell(length(pt(whichPt).sz),1);


for j = 1:length(pt(whichPt).sz)
    
    %% Define time windows
    totalTime = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j).runTimes(1,1);
    nchunks = ceil(totalTime/window);
    
    gdf_all = [];
    
    %% Get all spikes in one big array for each seizure
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
            gdf = removeEKGArtifact(gdf,gdf_ekg,prox);
        end
        
        gdf_all = [gdf_all;gdf];
        
        if rmDepth == 1
            gdf_all = removeChs(gdf_all,electrodeData,rmType);
        end
        
    end
    
    chunk_spikes = zeros(nchunks,1);
    chunk_spikes_chs = zeros(nchunks,nchs);
    
    %% Populate spike frequency for the time chunk
    for tt =  1:nchunks
        
        times = [(tt-1)*window + pt(whichPt).sz(j).runTimes(1,1),...
            tt*window+pt(whichPt).sz(j).runTimes(1,1)];
        
         % get the appropriate spikes in this time
         correct_spikes =  gdf(gdf(:,2) >= times(1) && gdf(:,2) <= times(2),:);
         
         % put the number of spikes into the chunk_spikes array
         chunk_spikes(tt) = size(correct_spikes,1);
         
         % get the number of spikes per channel
         for k = 1:size(correct_spikes,1)
            ch = correct_spikes(k,1);
            chunk_spikes_chs(tt,k) = chunk_spikes_chs(tt,k) + 1;
             
         end
        
    end
    
    % Divide by the size of the window to get the frequency
    chunk_spikes = chunk_spikes/window;
    chunk_spikes_chs = chunk_spikes_chs/window;
    
    % Populate the cell array containing spike freq for all seizures
    spike_freq{j} = chunk_spikes;
    spike_freq_ch{j} = chunk_spikes_chs;
      
end

%% Now find a way to combine the different seizures for easy plotting
all_spike_freq{1} = spike_freq{1};
all_spike_freq_ch{1} = spike_freq_ch{1};

for j = 2:length(pt(whichPt).sz)
   szTime = pt(whichPt).sz(j).onset;
   szTimeLast = pt(whichPt).sz(j-1).onset;
   
   % if the last seizure time is within 12 hours of the current seizure
   % time, append the appropriate times and look at it all as one
   % continuous chunk
   if szTime - szTimeLast <= totalTime/2
       
       % We only need to get the times from the end of the last run to the
       % end of this run
       timeToAdd = pt(whichPt).sz(j).runTimes(end,2) - pt(whichPt).sz(j-1).runTimes(end,2);
       
       all_spike_freq{j-1} = [all_spike_freq{j-1};***];
       
   end
   
    
end


end