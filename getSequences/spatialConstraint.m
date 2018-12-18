%{
 
   Perform spatial constraint: require that each spike in the sequence be
   close enough to the previous spike.

   There is an additional caveat that if the spikes are far apart, but it
   is very common that Channel 1 has a spike right before Channel 2, then I
   allow it anyway.

%}

function part = spatialConstraint(temp, xyChan, ...
    minConcurrentFreqAbs, minConcurrentFreqRel, maxSpeed, fs, minSpikesCloseEnough,doMorph)

part = [];
n_chans = size(xyChan,1);

%% Make Network of pairwise co-occurances
network = zeros(n_chans,n_chans);

if doMorph == 1
    jump = 4;
else
    jump = 2;
end

% Loop over sequences
for c = 1:jump:size(temp,2)
    
    % Get column of channels in a specific sequence
    col   = temp(:,c); 
    
    % remove zeros
    col(col==0) = [];  
    
    % Loop through all channels that have a spike in that sequence
    for r = 2:length(col)        
        
        % add a 1 to this network matrix if there is an occasion where
        % channel1 immediately preceded channel2
        network(col(r-1),col(r)) = network(col(r-1),col(r)) + 1;
    end  
end

%% Create vector of all nonzero entries in temp matrix (used later)


% This is just an array of channel numbers corresponding to the spikes for
% each spike sequence, starting with the next channel after the first. So
% it is nxm, where n is the number of spikes in the sequence (padded with
% zeros) and m is the number of sequences.
entries  = temp(2:end,1:jump:end);


% This reshapes entries from an nxm (number of spikes x number
% of sequences) array to a 1x(n*m) array, where it goes Channel according
% to first spike from first sequence, then channel from first spike from
% 2nd sequence, ..., then channel from 2nd spike from 1st sequence, then
% channel from 2nd spike from 2nd seq, etc. I think the point is just to
% get a single dimension array to very easily count (down below) how often
% a given channel shows up in the head channels sequences
entries = reshape(entries.',1,numel(entries));

% remove zeros
entries(entries==0)=[];

% for testing
%speedForTesting = [];

%% Test whether each spike falls in legal location
for c = 1:jump:size(temp,2)
    % loop through the sequences
    
    
    col             = temp(:,c); % which channel
    ticks           = temp(:,c+1); % which time
    if doMorph == 1
        height = temp(:,c+2);
        width = temp(:,c+3);
    end
    col(col==0)     = [];
    ticks(ticks==0) = [];
    if doMorph == 1
        height(height==0)     = [];
        width(width==0)     = [];
    end
    r               = 2; % the row (which spike in the sequence). Start with the second spike
    
    % loop through the spikes in the sequence
    while r <= size(col,1)
        
        % For the previous spike, get the x, y, z location of the channel
        prevChan = col(r-1,1);
        xyzPrevChan = xyChan(prevChan,2:end);
        
        %Get the time of the previous spike
        timePrevSpike = ticks(r-1,1);
        
        % Get the x,y,z location of the current channel
        currChan = col(r,1);
        xyzCurrChan = xyChan(currChan,2:end);
        
        % Get the time of the current spike
        timeCurrSpike = ticks(r,1);
        
        % Calculate the necessary speed the spike needed to travel to go
        % from one electrode to the other, and whether that's allowable
        distTravel = sqrt(sum((xyzPrevChan-xyzCurrChan).^2));
        timeTravel = timeCurrSpike-timePrevSpike;
        
        if timeTravel == 0
            % if the spike occured at the same time, assume the travel time
            % was actually one over the sampling rate (so the shortest time
            % step)
           timeTravel =  1/fs;
        end
        
        speedTravel = distTravel/timeTravel;
        %speedForTesting = [speedForTesting,speedTravel];
        
        % if the 2 channels are NOT within the allowable distance from each
        % other
        if speedTravel > maxSpeed
            

            % Test for frequency of connection
            freq = network(col(r-1,1),col(r,1));
            
            % If this frequency is too low (either in absolutes or relative
            % to how connected the prior channel is in general
            if freq < minConcurrentFreqAbs ||...
                    freq/sum(network(col(r-1,1),:)) < minConcurrentFreqRel
                
                % Remove spike from sequence. Note that I have not yet
                % removed the subsequent spikes in the sequence
                col(r)   = [];
                ticks(r) = [];
                if doMorph == 1
                    height(r) = [];
                    width(r) = [];
                end
                    
                r = r+0;        % do not increment row
                
            else
                
                r = r+1;        % move to next row
                
            end
            
        else
            
            r = r + 1;          % move to next row
            
        end
    end
    
    
    % If, after this, there are still enough spikes in the sequence,
    % include the sequence
    if length(col) >= minSpikesCloseEnough
        % Add dummy zeros, concatenate
        diff = length(col) - size(part,1);
        if diff < 0
            col = vertcat(col, zeros(abs(diff),1));
            ticks = vertcat(ticks, zeros(abs(diff),1));
            if doMorph == 1
                height = vertcat(height, zeros(abs(diff),1));
                width = vertcat(width, zeros(abs(diff),1));
            end
        elseif diff > 0
            part = vertcat(part, zeros(diff,size(part,2)));
        end
        if doMorph == 1
            part = [part, col, ticks,height,width];
        else
            part = [part, col, ticks];
        
        end
    end
end

end