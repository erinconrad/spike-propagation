function MI = getMI(avgRecruitmentLat,wij)

nChannels = length(avgRecruitmentLat);

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
    
    MI = outsidemultiple*numerator/denominator;
    



end