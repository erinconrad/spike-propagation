function Patient = getSpatialOrg(Patient,indexToms,dmin)


nSegments = size(Patient.sequences,2); % Number of segments
nChannels = size(Patient.xyChan,1); % Number of channels

%Initialize spatial org, just for each segment
spatialOrg = nan(1,nSegments);
avgRecruitmentLat = cell(1,nSegments);
allChannelAvgRecruitmentLat = nan(1,nSegments);
for iSeg = 1:nSegments
    avgRecruitmentLat{iSeg} = nan(1,nChannels);
    
end

% Calculate mean recruitment latency per channel
for iSeg = 1:nSegments
    for iChannel = 1:nChannels
        
        % Average across the rows of the recruitment latencies
        avgRows = mean(Patient.recruitmentLatencySingle{iSeg}{iChannel},1);
        
        % The first value of this is the average recruitment time (the
        % second is the average of the lead channel numbers which is
        % meaningless)
        avgRecruitmentLat{iSeg}(iChannel) = avgRows(1);
    end
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
           dij(iChannel,jChannel) =  sqrt(sum((Patient.xyChan(iChannel,2:end) - Patient.xyChan(jChannel,2:end)).^2));
           if dij(iChannel,jChannel) < dmin
               wij(iChannel,jChannel) = 1/dij(iChannel,jChannel);
               if wij(iChannel,jChannel) == inf
                   wij(iChannel,jChannel) = 0;
                   fprintf('Warning, channels %d and %d have the same location, making spatial weight 0\n',iChannel,jChannel);
               end
           else
               wij(iChannel,jChannel) = 0;
           end
       end
      
   end
    
end


% Calculate average mean recruitment latency across all channels
for iSeg = 1:nSegments
   allChannelAvgRecruitmentLat(iSeg) = nanmean(avgRecruitmentLat{iSeg});  
end

% Calculate spatial organization
for iSeg = 1:nSegments
    outsidemultiple = nChannels/nansum(nansum(wij));
    numerator = 0;
    denominator = 0;
    
    
    % calculate numerator
    for iChannel = 1:nChannels
        if isnan(avgRecruitmentLat{iSeg}(iChannel)) == 1 
            continue
        end
        for jChannel = 1:nChannels
            if jChannel == iChannel
                continue
            end
            if isnan(avgRecruitmentLat{iSeg}(jChannel)) == 1 
                continue
            end
            numerator = numerator + wij(iChannel,jChannel)*...
                (avgRecruitmentLat{iSeg}(iChannel)-allChannelAvgRecruitmentLat(iSeg))*...
                (avgRecruitmentLat{iSeg}(jChannel)-allChannelAvgRecruitmentLat(iSeg));
        end
    end
    
    % calculate denominator
    for iChannel = 1:nChannels
        if isnan(avgRecruitmentLat{iSeg}(iChannel)) == 1 
            continue
        end
        denominator = denominator + (avgRecruitmentLat{iSeg}(iChannel)-allChannelAvgRecruitmentLat(iSeg))^2;       
    end
    
    spatialOrg(iSeg) = outsidemultiple*numerator/denominator;
    
end

Patient.spatialOrg = spatialOrg;
Patient.avgRecruitmentLat = avgRecruitmentLat;

end