function power = alphaDelta(whichPts)

%% Load file paths, etc.
[~,~,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];

seq_file = 'long_seq.mat';
power_file = 'power.mat';

load([structFolder,seq_file])

if exist([structFolder,power_file],'file') ~= 0
    load([structFolder,power_file])
    
else
    power = struct;
end


if isempty(whichPts) == 1
    for i = 1:length(pt)
        if isempty(pt(i).seq_matrix) == 0
            whichPts = [whichPts,i];
        end
    end
end

for whichPt = whichPts
    fprintf('Doing %s\n',pt(whichPt).name)
    if size(power,1) >= whichPt
        if isfield(power(whichPt),'ad_rat_fft') == 1
            if isempty(power(whichPt).ad_rat_fft) == 0
                fprintf('Already did %s, skipping\n',pt(whichPt).name);
                continue
            end
        end
    else
    
    end
    
    fs = pt(whichPt).fs;
    dataName = pt(whichPt).ieeg_name;
    channels = pt(whichPt).channels;
    nch = length(channels);
    
    nbands = ceil(fs/2/10);
    
    szTimes = pt(whichPt).newSzTimes;
    
    
    ad_rat = zeros(nch,size(pt(whichPt).runTimes,1));
    ad_rat_band = zeros(nch,size(pt(whichPt).runTimes,1));
    all_p = zeros(nch,nbands,size(pt(whichPt).runTimes,1));
    times_out = mean(pt(whichPt).runTimes,1);
    
    %  Loop over run times
    for tt = 1:size(pt(whichPt).runTimes,1)
        fprintf('Doing chunk %d of %d\n',tt,size(pt(whichPt).runTimes,1));
        
        % Get the desired indices
        desiredTimes = pt(whichPt).runTimes(tt,:);
        
        % get times to clip
        clipTime = [-1*60 0];
        szTimesPlusClip = szTimes + repmat(clipTime,size(szTimes,1),1);
        szTimesT = szTimesPlusClip';
        out=range_intersection(desiredTimes,szTimesT(:));
        if isempty(out) == 1
            indicesToClip = []; 
        else
            indicesToClip = round(out(1)*fs):round(out(2)*fs);
        end
        
        %desiredTimes = desiredTimesFake(tt,:);
        indices = round(desiredTimes(1)*fs):round(desiredTimes(2)*fs);
        
        
        
        % Get the data
        
        data = getiEEGData(dataName,channels,indices,pwfile);
        
        fprintf('Retrieved data, doing analysis.\n');
        
        
        % remove nans
        data.values(isnan(data.values)) = 0;
        
        % Remove seizure times
        data.values(indicesToClip,:) = [];
        
        
        % Loop over channels
        
        for dd = 1:nch
            X = data.values(:,dd);
            
            % subtract mean
            X = X - mean(X);
            
            
            %% fft approach
           
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
            
            %fprintf('Got FFT approach\n');
            
            % Get all frequency bands
            %{
            for bb = 1:nbands
                frange = (bb-1)*10:min(bb*10,fs/2);
                all_p(dd,bb,tt) = sum(P(freqs>=frange(1) & freqs<=frange(2)));
                
            end
            %}

            %% Bandpass approach
           
         
            
            ad_rat_band(dd,tt) = bandpower(X,fs,[8 13])/bandpower(X,fs,[1 4]);
           
            
            
            
        end
        
        %fprintf('Finished analysis\n');
        
        
       
    
    end
    
  %  power(whichPt).all_p = all_p;
    power(whichPt).ad_rat_fft = ad_rat;
    power(whichPt).ad_rat_band = ad_rat_band;
    power(whichPt).times = times_out;
    save([structFolder,power_file])
end


end