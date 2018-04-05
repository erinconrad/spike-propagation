%{

This function takes information about recruitment latency and calculates
the spatial organization metric

%}

function [avgRecruitmentLat,spatialOrg] = getSpatialOrg(recruitmentLatencySingle,xyChan,dmin)


nChannels = size(xyChan,1); % Number of channels

%Initialize average recruitment latency
avgRecruitmentLat = nan(1,nChannels);


% Calculate mean recruitment latency per channel
for iChannel = 1:nChannels

    % Average across the rows of the recruitment latencies
    avgRows = mean(recruitmentLatencySingle{iChannel},1);

    % The first value of this is the average recruitment time (the
    % second is the average of the lead channel numbers which is
    % meaningless)
    avgRecruitmentLat(iChannel) = avgRows(1);
end


% Calculate dij, distances between channels and weights wij
dij = nan(nChannels,nChannels);
wij = nan(nChannels,nChannels);
for iChannel = 1:nChannels
   for jChannel = 1:nChannels
       if iChannel == jChannel
           dij(iChannel,jChannel) = nan;
           wij(iChannel,jChannel) = nan;
       else
           % calculate Euclidean distance between the two channels
           dij(iChannel,jChannel) =  sqrt(sum((xyChan(iChannel,2:end) - xyChan(jChannel,2:end)).^2));
           
           % if this distance is less than a minimum distance
           if dij(iChannel,jChannel) <= dmin
               
               % make the weight be 1/d
               wij(iChannel,jChannel) = 1/dij(iChannel,jChannel);
               
               if wij(iChannel,jChannel) == inf
                   wij(iChannel,jChannel) = 0;
                   fprintf('Warning, channels %d and %d have the same location, making spatial weight 0\n',iChannel,jChannel);
               end
           else
               
               % if they're further, set the weight to 0
               wij(iChannel,jChannel) = 0;
           end
       end
      
   end
    
end


% Calculate average mean recruitment latency across all channels
allChannelAvgRecruitmentLat = nanmean(avgRecruitmentLat);  


% Calculate spatial organization
    outsidemultiple = nChannels/nansum(nansum(wij));
    numerator = 0;
    denominator = 0;
    
    
    % calculate numerator
    for iChannel = 1:nChannels
        if isnan(avgRecruitmentLat(iChannel)) == 1 
            continue
        end
        for jChannel = 1:nChannels
            if jChannel == iChannel
                continue
            end
            if isnan(avgRecruitmentLat(jChannel)) == 1 
                continue
            end
            numerator = numerator + wij(iChannel,jChannel)*...
                (avgRecruitmentLat(iChannel)-allChannelAvgRecruitmentLat)*...
                (avgRecruitmentLat(jChannel)-allChannelAvgRecruitmentLat);
        end
    end
    
    % calculate denominator
    for iChannel = 1:nChannels
        if isnan(avgRecruitmentLat(iChannel)) == 1 
            continue
        end
        denominator = denominator + (avgRecruitmentLat(iChannel)-allChannelAvgRecruitmentLat)^2;       
    end
    
    spatialOrg = outsidemultiple*numerator/denominator;
    

end