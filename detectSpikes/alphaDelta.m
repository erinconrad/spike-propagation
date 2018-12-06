function power = alphaDelta(pt,whichPts)

%% Load file paths, etc.
[~,~,~,~,pwfile] = fileLocations;
power = struct;
%t = 323420.99;
%desiredTimesFake = [323420.99 323420.99+15;323420.99+16 323420.99+30];

for whichPt = whichPts
    
    fs = pt(whichPt).fs;
    dataName = pt(whichPt).ieeg_name;
    channels = pt(whichPt).channels;
    nch = length(channels);
    
    nbands = ceil(fs/2/10);
    
    
    ad_rat = zeros(nch,size(pt(whichPt).runTimes,1));
    all_p = zeros(nch,nbands,size(pt(whichPt).runTimes,1));
    times_out = mean(pt(whichPt).runTimes,1);
    
    %  Loop over run times
    for tt = 1:size(pt(whichPt).runTimes,1)
        fprintf('Doing chunk %d of %d\n',tt,size(pt(whichPt).runTimes,1));
        
        % Get the desired indices
        desiredTimes = pt(whichPt).runTimes(tt,:);
        
        %desiredTimes = desiredTimesFake(tt,:);
        indices = round(desiredTimes(1)*fs):round(desiredTimes(2)*fs);
        
        % Get the data
        tic
        data = getiEEGData(dataName,channels,indices,pwfile);
        toc
        fprintf('Retrieved data\n');
        
        % remove nans
        data.values(isnan(data.values)) = 0;
        
        % Loop over channels
        
        for dd = 1:nch
            X = data.values(:,dd);
            
            % subtract mean
            X = X - mean(X);
            
            
            
            % Calculate fft
            Y = fft(X);
            
            % Get power
            P = abs(Y).^2;
            freqs = linspace(0,fs,length(P)+1);
            freqs = freqs(1:end-1);
            
            % Take first half
            P = P(1:ceil(length(P)/2));
            freqs = freqs(1:ceil(length(freqs)/2));
            
            %plot(freqs,P);
            
            % Get alpha/delta ratio
            alpha = sum(P(freqs>=8 & freqs<=13));
            delta = sum(P(freqs>=1 & freqs<=4));
            ad_rat(dd,tt) = alpha/delta;
            
            % Get all frequency bands
            for bb = 1:nbands
                frange = (bb-1)*10:min(bb*10,fs/2);
                all_p(dd,bb,tt) = sum(P(freqs>=frange(1) & freqs<=frange(2)));
                
            end
            
        end
        
        
       
    
    end
    
    power(whichPt).all_p = all_p;
    power(whichPt).ad_rat = ad_rat;
    power(whichPt).times = times_out;
end


end