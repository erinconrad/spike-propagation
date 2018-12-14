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
    if size(power,1) >= whichPt
        if isfield(power(whichPt),'ad_rat_fft') == 1
            if isempty(power(whichPt).ad_rat_fft) == 0
                fprintf('Already did %s, skipping\n',pt(whichPt).name);
                continue
            end
        end
    else
    fprintf('Doing %s\n',pt(whichPt).name)
    end
    
    fs = pt(whichPt).fs;
    dataName = pt(whichPt).ieeg_name;
    channels = pt(whichPt).channels;
    nch = length(channels);
    
    nbands = ceil(fs/2/10);
    
    
    ad_rat = zeros(nch,size(pt(whichPt).runTimes,1));
    ad_rat_band = zeros(nch,size(pt(whichPt).runTimes,1));
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
            
            % Get all frequency bands
            for bb = 1:nbands
                frange = (bb-1)*10:min(bb*10,fs/2);
                all_p(dd,bb,tt) = sum(P(freqs>=frange(1) & freqs<=frange(2)));
                
            end
            

            %% Bandpass approach
            alpha_freq = bandpass(X,[8 13],fs);
            alpha_pow = mean(alpha_freq.^2);

            delta_freq = bandpass(X,[1 4],fs);
            delta_pow = mean(delta_freq.^2);
            ad_rat_band(dd,tt) = alpha_pow./delta_pow;

        end
        
        
       
    
    end
    
  %  power(whichPt).all_p = all_p;
    power(whichPt).ad_rat_fft = ad_rat;
    power(whichPt).ad_rat_band = ad_rat_band;
    power(whichPt).times = times_out;
    save([structFolder,power_file])
end


end