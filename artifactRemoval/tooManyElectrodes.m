%{
This function loops through the gdf array and discards spikes that occur in
too high a percentage of channels within too close a time to each other

%}


function gdf = tooManyElectrodes(gdf,maxChannels,multiChTime)

discard = zeros(size(gdf,1),1);
i = 1;
% Loop through the spikes (we are going to variably move through the spikes)
while 1
    scount = 1;

    % if i is the last spike, break
    if i >= size(gdf,1)
        break
    end
    
    % the time we are comparing other times to
    refTime = gdf(i,2);

   % for each spike, loop through subsequent spikes to see how many
   % there are across multiple channels occuring at the same time
   for j = i+1:size(gdf,1)

       % Get the time difference between the first spike and the last
       % spike
       tdiff = gdf(j,2)-refTime;

       % If the difference is small enough, increase the count
       if tdiff < multiChTime
           scount = scount + 1;
       else
           % break out of the inner loop and go to the next spike to start counting 
           break
       end

       % if the total number of channels spiking in this very close
       % proximity is >max allowed
       if scount > maxChannels

            % Mark them to be discarded
            discard(i:j) = ones(length(j-i+1),1);

            % Break the inner for loop and move to the next spike
            i = j;
            break
       end

   end

   % advance the index of the spike
   i = i+1;
end

gdf(find(discard==1),:) = [];

end