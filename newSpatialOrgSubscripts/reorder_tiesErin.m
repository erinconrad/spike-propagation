%{
 
    Function: detect sequences in gdf. Look for ties- reorder tied channels
    by proximity to closest untied event. 
          
%}
 
function [overall] = reorder_tiesErin(overall, xyChan)

for c = 2:2:size(overall,2)     % loop ticks columns in overall
    
    ticksCol = overall(:,c);
    chanCol  = overall(:,c-1);
    
    ticksCol(ticksCol==0)=[]; % remove zeros resulting from concactenation
    chanCol(chanCol==0)=[];   
    
    u = unique(ticksCol)';
    for tick = u
        ind  = find(tick == ticksCol);
        if length(ind) > 1    % this means there are ties for this tick
            
            % if not the first step, take proximity to prev entry
            if ind(1) ~= 1
                
                prevInx        = ind(1)-1;
                prevChan       = chanCol(prevInx,1);
                prevX = xyChan(prevChan,2);
                prevY = xyChan(prevChan,3);
                prevZ = xyChan(prevChan,4);
                delta        = [];
                
                for currInx  = ind
                    currChan = chanCol(currInx,1);
                    currX = xyChan(currChan,2);
                    currY = xyChan(currChan,3);
                    currZ = xyChan(currChan,4);
                    dist     = sqrt(((prevX-currX).^2)+((prevY-currY).^2)+(prevZ-currZ).^2);
                    delta    = [delta, dist];
                end
                
                % Reorder, update overall
                [Y,I] = sort(delta);
                I = I + prevInx;    % Offset 
                chanCol(ind,1) = chanCol(I,1);
                
                % Add the dummy zeros back
                z = size(overall,1)-size(chanCol,1);
                chanCol = vertcat(chanCol, zeros(z,1));
                ticksCol = vertcat(ticksCol, zeros(z,1));
                overall(:,c)   = ticksCol;
                overall(:,c-1) = chanCol;
                
            % If tie is first step- go by proximity to first untied entry
            elseif ind(end) < size(chanCol,1)
                
                postInx  = ind(end)+1;
                postChan = chanCol(postInx,1);
                postX = xyChan(postChan,2);
                postY = xyChan(postChan,3);
                postZ = xyChan(postChan,4);
                delta    = [];
                
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
                I = flipud(I);  % arrange from furthest to closest
                chanCol(ind,1) = chanCol(I,1);
                
                % Add the dummy zeros back
                z = size(overall,1)-size(chanCol,1);
                chanCol = vertcat(chanCol, zeros(z,1));
                ticksCol = vertcat(ticksCol, zeros(z,1));
                overall(:,c)   = ticksCol;
                overall(:,c-1) = chanCol;
                
            end
        end
    end
    
end

end