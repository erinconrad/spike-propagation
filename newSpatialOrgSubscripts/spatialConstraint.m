%{
 
   Perform spatial constraint: require that each spike in the sequence be
   close enough to the previous spike

%}

function part = spatialConstraint(temp, xyChan)

minConcurrentFreq = 0.15;
maxDist = 4;

part = [];

%% Create vector of all nonzero entries in temp matrix (used later)

% each row is a spike (or zero padding) and then each column is a channel
% number going through all the different sequences
entries  = temp(2:end,1:2:end);

% not sure what this does
entries = reshape(entries.',1,numel(entries));

entries(entries==0)=[];

%% Test whether each spike falls in legal location
for c = 1:2:size(temp,2)
    % loop through the sequences
    
    
    col             = temp(:,c); % which channel
    ticks           = temp(:,c+1); % which time
    col(col==0)     = [];
    ticks(ticks==0) = [];
    r               = 2; % the row (which spike in the sequence). Start with the second spike
    
    % loop through the spikes in the sequence
    while r <= size(col,1)
        
        % For the previous spike, get the x, y, z location of the channel
        prevChan = col(r-1,1);
        xyzPrevChan = xyChan(prevChan,2:end);
        
        % Get the x,y,z location of the current channel
        currChan = col(r,1);
        xyzCurrChan = xyChan(currChan,2:end);
        
        % if the 2 channels are NOT within the allowable distance from each
        % other
        if sqrt((xyzPrevChan-xyzCurrChan).^2) < maxDist

            % Test for frequency in temp matrix. This is a test of how
            % often the current test channel shows up in the spike trains
            % for the current "head" channel defining the temp matrix
            occurances = find(col(r,1)==entries);
            
            % If this channel doesn't show up enough in the head channel's
            % spike trains
            if length(occurances) / length(entries) < minConcurrentFreq
                
                % Remove spike from sequence
                col(r)   = [];
                ticks(r) = [];
                r = r+0;        % do not increment row
                
            else
                
                r = r+1;        % move to next row
                
            end
            
        else
            
            r = r + 1;          % move to next row
            
        end
    end
    
    
    % If there are more than 2 spikes??
    if length(col) >= 3
        % Add dummy zeros, concactenate
        diff = length(col) - size(part,1);
        if diff < 0
            col = vertcat(col, zeros(abs(diff),1));
            ticks = vertcat(ticks, zeros(abs(diff),1));
        elseif diff > 0
            part = vertcat(part, zeros(diff,size(part,2)));
        end
        
        part = [part, col, ticks];
    end
end

end