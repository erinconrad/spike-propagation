function ad_remove_spikes(whichPts)

%% Parameters
alpha_freq = [8 13];
delta_freq = [1 4];
sp_surround = [-0.5 0.5];

%% Load file paths, etc.
[~,~,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];

seq_file = 'long_seq.mat';
power_file = 'power_check.mat';

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
    fs = pt(whichPt).fs;
    dataName = pt(whichPt).ieeg_name;
    channels = pt(whichPt).channels;
    nch = length(channels);
    szTimes = pt(whichPt).newSzTimes;
    
    power(whichPt).times = mean(pt(whichPt).runTimes,2);
    if isfield(power(whichPt),'alpha') == 0 || isempty(power(whichPt).alpha) == 1
        power(whichPt).alpha.rm_spike = zeros(nch,size(pt(whichPt).runTimes,1));
        power(whichPt).delta.rm_spike = zeros(nch,size(pt(whichPt).runTimes,1));
        
        %power(whichPt).ad_rat = zeros(nch,size(pt(whichPt).runTimes,1));
        %power(whichPt).ad_rat_alt = zeros(nch,size(pt(whichPt).runTimes,1));
        
        power(whichPt).alpha.keep_spike = zeros(nch,size(pt(whichPt).runTimes,1));
        power(whichPt).delta.keep_spike = zeros(nch,size(pt(whichPt).runTimes,1));
        
        power(whichPt).finished = zeros(size(pt(whichPt).runTimes,1),1);
    end
    
    % Get the spike times (I am taking all times in the seq_matrix)
    seq_matrix = pt(whichPt).seq_matrix;
    all_spike_times = seq_matrix(~isnan(seq_matrix));
    all_spike_chs = repmat([1:size(seq_matrix,1)],size(seq_matrix,2),1)';
    all_spike_chs = all_spike_chs(~isnan(seq_matrix));
    all_spikes = [all_spike_chs,all_spike_times];
    
    
    
    %  Loop over run times
    for tt = 1:size(pt(whichPt).runTimes,1)
        
        if power(whichPt).finished(tt) == 1
            fprintf('Already did chunk %d for %s, skipping\n',...
                tt,pt(whichPt).name);
            continue
        end
        
        fprintf('Doing chunk %d of %d for %s\n',tt,...
            size(pt(whichPt).runTimes,1),pt(whichPt).name);
        
        
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
        indices = max(round(desiredTimes(1)*fs),1):round(desiredTimes(2)*fs);
        
        % get the data
        data = getiEEGData(dataName,channels,indices,pwfile);
        
        % remove nans
        data.values(isnan(data.values)) = 0;

        % Remove seizure times
        data.values(max(indicesToClip-indices(1),1),:) = [];
        nch = length(channels);
        
        %% Get the indices corresponding to the spikes
        
        % find the spike times that are in the current interval
        curr_spikes = all_spikes(all_spikes(:,2) > desiredTimes(1) & all_spikes(:,2) < desiredTimes(2),:);
        
        % convert to indices
        curr_spikes(:,2) = (curr_spikes(:,2) - desiredTimes(1))*fs;
        
        
        for dd = 1:nch
            
            % get values for that channel
            X = data.values(:,dd);
            X = X - mean(X);
            
            % Find the spikes in that channel
            spikes_in_ch = curr_spikes(curr_spikes(:,1) == dd,:);
            
            % Plot the spikes as an example
            if 0
            for i = 1:size(spikes_in_ch,1)
                figure
                times = round(spikes_in_ch(i,2) - 7*fs):round(spikes_in_ch(i,2) + 7*fs);
                plot(times,X(times))
                hold on
                plot(round(spikes_in_ch(i,2)),X(round(spikes_in_ch(i,2))),'ro')
                pause
                close(gcf)
            end
            end
            
            % Bandpass filter the data
            alpha = bandpass(X,alpha_freq,fs,'ImpulseResponse','iir');
            delta = bandpass(X,delta_freq,fs,'ImpulseResponse','iir');
            
            % Get the power
            alpha_power = alpha.^2;
            delta_power = delta.^2;
            
            % Get spike exclusion indices
            sp_exclusion_idx = zeros(size(spikes_in_ch,1),sp_surround(2)*fs*2+1);
            for i = 1:size(spikes_in_ch,1)
                sp_exclusion_idx(i,:) = round(spikes_in_ch(i,2) + ...
                    sp_surround(1)*fs): round(spikes_in_ch(i,2) + ...
                    sp_surround(2)*fs);
            end

            
            % Look at the signal and plot the spike exclusion times
            if 0
                figure
                Y = delta;
                plot(Y,'k')
                hold on
                for i = 1:size(sp_exclusion_idx,1)
                    plot(sp_exclusion_idx(i,:),Y(sp_exclusion_idx(i,:)),'r');
                end
                pause
                close(gcf)
            end
            
            % Sum the power
            keep_idx = ones(length(alpha_power),1);
            keep_idx(sp_exclusion_idx) = 0;
            
            alpha_power_sum_ex = sum(alpha_power(keep_idx));
            alpha_power_sum_in = sum(alpha_power);
            
            delta_power_sum_in = sum(delta_power);
            delta_power_sum_ex = sum(delta_power(keep_idx));
            
            power(whichPt).alpha.rm_spike(dd,tt) = alpha_power_sum_ex;
            power(whichPt).alpha.keep_spike(dd,tt) = alpha_power_sum_in;
            
            power(whichPt).delta.rm_spike(dd,tt) = delta_power_sum_ex;
            power(whichPt).delta.keep_spike(dd,tt) = delta_power_sum_in;
            
        end
        
        power(whichPt).finished(tt) = 1;
        save([structFolder,power_file],'power')
        
    end
    
end


end