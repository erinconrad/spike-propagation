%{
 
    Function: Loop through sequences in gdf. Look for ties and reorder tied channels
    by proximity to closest untied event. 
          
%}
 
function [overall] = reorder_tiesErin(overall, xyChan,doMorph)

if doMorph == 1, jump = 4; else, jump = 2; end

% Loop through the even columns in overall, which contain the times of the
% spikes for each of the sequences. We are essentially looping over
% sequences
for c = 2:jump:size(overall,2)     
    
     % the column with the times for each spike in the current sequence
    ticksCol = overall(:,c);
    
    % the column with the channels for each spike in the current sequence
    chanCol  = overall(:,c-1);
    
    if doMorph == 1
        heightCol = overall(:,c+1);
        widthCol = overall(:,c+2);
    end
    
    % remove zeros resulting from concatenation
    ticksCol(ticksCol==0)=[]; 
    chanCol(chanCol==0)=[];   
    
    if doMorph == 1
        heightCol(ticksCol==0) = [];
        widthCol(ticksCol==0) = [];
    end
   
    
    % Find the unique times
    u = unique(ticksCol)';
    
    % Loop through the unique times
    for tick = u
        
        % This is the number of spikes with the same time
        ind  = find(tick == ticksCol);
        
        if length(ind) > 1    % this means there are ties for this tick
            
            % if not the first step, take proximity to prev entry. So this
            % is getting the location of the last non-tied channel, or the
            % last unique entry
            if ind(1) ~= 1
                
                prevInx        = ind(1)-1;
                prevChan       = chanCol(prevInx,1);
                prevX = xyChan(prevChan,2);
                prevY = xyChan(prevChan,3);
                prevZ = xyChan(prevChan,4);
                
                
                % This will be an array of distances between all tied
                % channels and the last non-tied channel
                delta        = [];
                
                % Loop through all the tied channels
                for currInx  = ind
                    currChan = chanCol(currInx,1);
                    currX = xyChan(currChan,2);
                    currY = xyChan(currChan,3);
                    currZ = xyChan(currChan,4);
                    % Location of current channel
                    
                    dist     = sqrt(((prevX-currX).^2)+((prevY-currY).^2)+(prevZ-currZ).^2);
                    delta    = [delta, dist];
                end
                
                % Sort the channels by proximity to the last non-tied
                % channel
                [Y,I] = sort(delta);
                I = I + prevInx;    % Offset 
                
                % Reorder the tied channels according to the sort
                chanCol(ind,1) = chanCol(I,1);
                
                if doMorph == 1
                    heightCol(ind,1) = heightCol(I,1);
                    widthCol(ind,1) = widthCol(I,1);
                end
                
            % If tie is first step- go by proximity to first untied entry.
            % This obviously only works if the whole thing isn't ties.
            elseif ind(end) < size(chanCol,1)
                
                % Now find the next untied entry
                postInx  = ind(end)+1;
                postChan = chanCol(postInx,1);
                postX = xyChan(postChan,2);
                postY = xyChan(postChan,3);
                postZ = xyChan(postChan,4);
                delta    = [];
                
                % Loop through the tied channels and make an array, delta,
                % that has the distances between the tied channels and the
                % next channel that doesn't tie.
                for currInx  = ind
                    currChan = chanCol(currInx,1);
                    currX = xyChan(currChan,2);
                    currY = xyChan(currChan,3);
                    currZ = xyChan(currChan,4);
                    dist     = sqrt(((postX-currX).^2)+((postY-currY).^2)+(postZ-currZ).^2);
                    delta    = [delta, dist];
                end
                
                % Reorder, update overall
                [Y,I] = sort(delta); % No offset needed b/c at start
                
                 % arrange from furthest to closest, because now we want
                 % the first tied channel to be far away from the next
                 % untied channel, and move progressively closer throughout
                 % the ties
                I = flipud(I); 
                
                % update the channel column
                chanCol(ind,1) = chanCol(I,1);
                
                if doMorph == 1
                    heightCol(ind,1) = heightCol(I,1);
                    widthCol(ind,1) = widthCol(I,1);
                end
               
                
            end
        end
    end
    
     % Add the dummy zeros back
    z = size(overall,1)-size(chanCol,1);
    chanCol = vertcat(chanCol, zeros(z,1));
    ticksCol = vertcat(ticksCol, zeros(z,1));
    
    if doMorph == 1
        z = size(overall,1)-size(heightCol,1);
        heightCol = vertcat(heightCol, zeros(z,1));
        z = size(overall,1)-size(widthCol,1);
        widthCol = vertcat(widthCol, zeros(z,1));
    end

    % update overall
    overall(:,c)   = ticksCol;
    overall(:,c-1) = chanCol;
    
    if doMorph ==1
        overall(:,c+1) = heightCol;
        overall(:,c+2) = widthCol;
    end
    
end

end