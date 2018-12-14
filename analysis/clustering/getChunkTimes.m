function [prop_pop,chunk_times,new_chunks] = getChunkTimes(allTimes,test_t,all_times_all,idx)

which_group = 1;
popular = mode(idx);
    
% fill up first one
new_chunks = [];
new_chunks{1} = [allTimes(1) allTimes(1) + test_t];

while 1
    curr_time(1) = new_chunks{end}(2);
    curr_time(2) = min(new_chunks{end}(2) + test_t,allTimes(which_group,2));
    new_chunks{end+1} = [curr_time(1) curr_time(2)];

    if curr_time(2) == allTimes(which_group,2)
        chunk_t = curr_time(2) - curr_time(1);
        which_group = which_group + 1;
        if which_group > size(allTimes,1) 
            break
        else
            new_chunks{end} = [new_chunks{end};...
                allTimes(which_group,1) allTimes(which_group,1) + test_t - chunk_t];
        end
    end
end

n_chunks = length(new_chunks);%ceil((max(all_times_all) - min(all_times_all))/test_t);
    prop_pop = zeros(n_chunks,1);
    chunk_times = zeros(n_chunks,1);
    
    for i = 1:n_chunks
        chunk_times(i) = new_chunks{i}(1,1);
        
        if size(new_chunks{i},1) == 1
            curr_times = new_chunks{i};
            % Get the spike indices in that time chunk
            chunk_spikes = find(all_times_all >= curr_times(1) & ...
                all_times_all <= curr_times(2));
            
        elseif size(new_chunks{i},1) == 2
            curr_times = new_chunks{i}(1,:);
            % Get the spike indices in that time chunk
            chunk_spikes = find(all_times_all >= curr_times(1) & ...
                all_times_all <= curr_times(2));
            
            curr_times = new_chunks{i}(2,:);
            % Get the spike indices in that time chunk
            chunk_spikes = [chunk_spikes;find(all_times_all >= curr_times(1) & ...
                all_times_all <= curr_times(2))];
            
        else
            error('wat');
        end
        
        %{
        
        % Get the time range for the chunk
        curr_times = [min(all_times_all) + (i-1)*test_t,...
           min(min(all_times_all) + i*test_t,max(all_times_all))];
       
        % Get the spike indices in that time chunk
        chunk_spikes = find(all_times_all >= curr_times(1) & ...
            all_times_all <= curr_times(2));
        %}
        
        % Get the proportion of spikes in the most popular cluster
        prop_pop(i) = sum(idx(chunk_spikes) == popular)/length(chunk_spikes);
       
    end


end