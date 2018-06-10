function [new_gdf,removed_for_noise,removed_as_empty] = removeNoisy(gdf_all,noise_all,empty_all)

% Remove spikes if they are this many seconds from "noisy" data
time_from_noisy = 10; %10 seconds

% I need the updated gdf files!!!!!!!!!

noise_thresh = nanmedian(noise_all(:,2))*10;

noisytimes = noise_all(noise_all(:,2)>noise_thresh,1);

noisy_idx = noise_all(:,2)>noise_thresh;
noise_diff = diff(noisy_idx);

noisy_chunks = [];
noisy_chunk_times = [];
switches = find(noise_diff~=0);

% This may fail if the end of the data is noisy!

% if the beginning of the data is not noisy
if noisy_idx(1) == 0
    % then the noisy periods will be the period starting with the first
    % switch to the 2nd switch, then the 3rd switch to the 4th switch, etc.
    for i = 1:2:length(switches)-1
        noisy_chunks = [noisy_chunks;switches(i) switches(i+1)];
        noisy_chunk_times = [noisy_chunk_times; noise_all(switches(i)) noise_all(switches(i+1))];
    end
% if the beginning of the data is noisy
elseif noisy_idx(1) == 1
    noisy_chunks = [1 switches(1)];
    for i = 2:2:length(switches)-1
       noisy_chunks = [noisy_chunks; switches(i) switches(i+1)];
       noisy_chunk_times = [noisy_chunk_times; noise_all(switches(i)) noise_all(switches(i+1))];
    end
end

rm_gdf = zeros(size(gdf_all,1),1);

%% Remove noisy data
for i = 1:size(gdf_all,1)
    
    for j = 1:size(noisy_chunk_times,1)
        
        
        % extend the time a little over which we will delete spikes
        temp_times(1) = max(noise_all(1,1),noisy_chunk_times(j,1) - time_from_noisy);
        temp_times(2) = min(noise_all(end,1),noisy_chunk_times(j,2) + time_from_noisy);
        
        % if the spike falls between those two times
        if gdf_all(i,2) >= temp_times(1) && gdf_all(i,2) <= temp_times(2)
            % mark it for deletion
            rm_gdf(i) = 1;
            break
        end
            
    end
    
end
removed_for_noise = find(rm_gdf == 1);
new_gdf = gdf_all;
% Delete spikes marked for deletion
new_gdf(rm_gdf == 1,:) = [];

rm_empty = zeros(size(new_gdf,1),1);

%% Remove empty data
for i = 1:size(new_gdf,1)
    for j = 1:size(empty_all,1)
       if abs(new_gdf(i,2)-empty_all) < 1
           rm_empty(i) = 1;
           break
       end
    end
    
end

removed_as_empty = find(rm_empty == 1);
new_gdf(rm_empty == 1,:) = [];

end