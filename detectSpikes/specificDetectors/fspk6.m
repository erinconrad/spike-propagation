%% Like fspk3 except I am changing parameters for depth electrodes

function [gdf,noise,removed] = fspk6(eeg,tmul,absthresh,n_chans,...
    srate,window,electrodes)

%% I CHANGED A BUNCH OF THINGS

%{
This program is the non-GUI version of the spike detection algorithm fspk.
It was broken down to run over the CHOP network remotely, not over Matlab.
USAGE:  out = fspk2('patient_num',tmul,absthresh)
INPUT:  'patient_num' <string>  -- Number of Patient to be analyzed
        tmul                       -- Threshold Multiplier
        absthresh                  -- Absolute Threshold
Outputs: .gdf files in same file as EEG -- array of [channel | tick ]
%}

% Check function input
if ~exist('tmul')
    fprintf('warning, no tmul entered, using 13\n');
    tmul=13;
end

if ~exist('absthresh')
    fprintf('warning, no absthresh entered, using 13\n');
    absthresh=300;
end

% Initialize parameters
rate   = srate;
chan   = 1:n_chans;

%% Depth electrodes are spikier
%fr     = 20;                 % high pass freq: probably don't have to change
%lfr    = 7;                  % low pass freq: probably don't have to change


spkdur = 220;                % spike duration must be less than this
spkdur = spkdur*rate/1000;   % convert to points;

%sharpdur = 200;  % duration of sharp wave
%sharpdur = sharpdur*rate/1000; % convert to points;


num_segs = ceil(size(eeg,1)/window);

% Initialize things
all_spikes  = [];
allout      = [];
overlap     = rate;
totalspikes = zeros(1, length(chan));
inc         = 0;
removed     = [];

% Read in eeg data
alldata     = eeg;            % timepnts x chans

% Check that it's enough time
if size(alldata,1) < window
    fprintf('Warning: should run spike detector on at least 60 seconds of data.\n');
end

% Initialize a noise matrix
noise = zeros(num_segs-1,n_chans);

% Iterate channels and detect spikes
for dd = 1:n_chans
    
    ch_type = electrodes(dd).type;
    
    fr     = 40;  % high pass freq, used to be 20
    lfr    = 7;   % low pass freq
    aftdur = 70;
    spikedur = 5;%5;
    fn_fr  = 10;
    
    %{
    if strcmp(ch_type,'D') == 1
        fr     = 30;%30;  % high pass freq, used to be 20
        lfr    = 7;   % low pass freq
        aftdur   = 70;%70;   % afterhyperpolarization wave must be longer than this
        spikedur = 10;%10;
        fn_fr = 1;

    else
        fr     = 20;  % high pass freq, used to be 20
        lfr    = 7;   % low pass freq
        aftdur = 150;
        spikedur = 20;
        fn_fr  = 1;
    end
    %}
    
    aftdur   = aftdur*rate/1000;   % convert to points;
    
    % Break the data into time segments, 1 minute each, so that the
    % threshold we are using to see if the spike rises above the background
    % is based on just the one minute we are considering.
    for tt = 1:num_segs
        time_points(1) = (tt-1)*window+1;
        time_points(2) = min(window*tt,size(eeg,1));
        
        if time_points(2) - time_points(1) < 100
            continue
        end
        
        if tt == num_segs
           time_points(2) = size(eeg,1);
        end
        
     
        
        out     = [];
        data    = alldata(time_points(1):time_points(2),dd);
        
  
        
        %% Skip spike detection if the amplitude is very low during the time period for that channel
        if sum(abs(data)) <= 1
            if isempty(removed) == 0 && removed(end,4) == 0 && removed(end,1) == dd
                removed(end,3) = time_points(2);
            else
                removed = [removed;dd time_points(1) time_points(2) 0];
            end
            continue
        end
        
        %% Run a noise detector; skip spike detection if it's a noisy minute for that channel
        % The second input is determining which method of noise detection
        % to use. If 1, then I use the sqrt of the sum of the squared
        % differences between adjacent time points. If 0, then I use the
        % RMS.
        noise_bin = findNoisyPeriods(data,1);
        noise(tt,dd) = noise_bin;

        if noise_bin == 1
            if isempty(removed) == 0 && removed(end,4) == 1 && removed(end,1) == dd
                removed(end,3) = time_points(2);
            else
                removed = [removed;dd time_points(1) time_points(2) 1];
            end
            continue
        end
        
        % Now filter
        %{
        f = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',srate);
           
        
        data = filtfilt(f,data);   
        %}
        
        
        
        %% re-adjust the mean of the data to be zero (if there is a weird dc shift)
        data = data - mean(data);
        
        %% Run the spike detector
        

      
        spikes   = [];
        shspikes = [];
        fndata   = eegfilt(data, fn_fr, 'hp',srate);

        % first look at the high frequency data for the 'spike' component
        HFdata    = eegfilt(fndata, fr, 'lp',srate);
       
        %{
        if strcmp(ch_type,'D') == 1
            HFdata = data;
        end
        %}
        
        lthresh = mean(abs(HFdata));  % this is the smallest the initial part of the spike can be
        thresh  = lthresh*tmul;     % this is the final threshold we want to impose
        sthresh = lthresh*tmul/3;   % this is the first run threshold

       
        [spp,spv] = FindPeaks(HFdata);

        idx      = find(diff(spp) <= spkdur);       % find the durations less than or equal to that of a spike
        startdx  = spp(idx);
        startdx1 = spp(idx+1);

        % check the amplitude of the waves of appropriate duration
        for i = 1:length(startdx)
            %spkmintic = spv(find(spv > startdx(i) & spv < startdx1(i)));  % find the valley that is between the two peaks
            spkmintic = spv((spv > startdx(i) & spv < startdx1(i)));
            %% commented out the second check
            if abs(HFdata(startdx1(i)) - HFdata(spkmintic)) > sthresh %& HFdata(startdx(i)) - HFdata(spkmintic) > lthresh  %#ok<AND2> % see if the peaks are big enough
                spikes(end+1,1) = spkmintic;                                  % add timestamp to the spike list
                spikes(end,2)   = (startdx1(i)-startdx(i))*1000/rate;         % add spike duration to list
                spikes(end,3)   = abs(HFdata(startdx1(i)) - HFdata(spkmintic));    % add spike amplitude to list
            end

        end

      
        % now have a list of spikes that have passed the 'spike' criterion.


        % try getting the sharp waves
        %idx = find(diff(spp) <= sharpdur & diff(spp) > spkdur);       % find the durations less than or equal to that of a spike
        %startdx = spp(idx);
        %startdx1 = spp(idx+1);

        % check the amplitude of the waves of appropriate duration
        %for i = 1:length(startdx)
        %    spkmintic = spv(find(spv > startdx(i) & spv < startdx1(i)));  % find the valley that is between the two peaks
        %    if HFdata(startdx1(i)) - HFdata(spkmintic) > sthresh & HFdata(startdx(i)) - HFdata(spkmintic) > lthresh  % see if the peaks are big enough
        %        shspikes(end+1,1) = spkmintic;                                  % if so add index to list
        %        shspikes(end,2) = HFdata(startdx1(i)) - HFdata(spkmintic);    % add spike amplitude to list
        %    end
        %end


        %spikes(:,1)        %these are the timestamps of the spikes
        %spikes(:,2)        %these are the durations in ms of the spike waves: at this point must already be less than spkdur
        %spikes(:,3)        %these are amplitudes in uV of the spike waves: at this point must already be = to or larger than sthresh
        spikes(:,4) = 0;    %these are the durations in ms of the afterhyperpolarization waves
        spikes(:,5) = 0;    %these are the amplitudes in uV of the afterhyperpolarization waves

        % now have a list of sharp waves that have passed criterion

        % check for after hyperpolarization
        dellist = [];
        
        

        LFdata = eegfilt(fndata, lfr, 'lp',srate);
        [hyperp,hyperv] = FindPeaks(LFdata);   % use to find the afterhyper wave
        olda = 0;  % this is for checking for repetitive spike markings for the same afterhyperpolarization
        for i = 1:size(spikes,1)
            % find the duration and amplitude of the slow waves, use this with the
            % amplitude of the spike waves to determine if it is a spike or not


            a = hyperp(find(hyperp > spikes(i,1)));          % find the times of the slow wave peaks following the spike

            try  % this try is just to catch waves that are on the edge of the data, where we try to look past the edge
                if a(2)-a(1) < aftdur                        % too short duration, not a spike, delete these from the list
                    dellist(end+1) = i;
                else 
                    % might be a spike so get the amplitude of the slow wave
                    spikes(i,4) = (a(2)-a(1))*1000/rate;       % add duration of afhp to the list
                    b = hyperv(find(hyperv > a(1) & hyperv < a(2))); % this is the valley
                    spikes(i,5) = abs(LFdata(a(1)) - LFdata(b));  % this is the amplitude of the afhp
                    if a(1) == olda    
                        % if this has the same afterhyperpolarization peak as the prev
                        %if strcmp(ch_type,'D') == 0
                            dellist(end+1) = i-1;           % spike then the prev spike should be deleted
                        %end
                    end
                end
                olda = a(1);

            catch
                dellist(end+1) = i;  % spike too close to the edge of the data
            end


        end

        s = spikes;
        
        
        
        spikes(dellist,:) = [];
        
       

        tooshort = [];
        toosmall = [];
        toosharp = [];

        % now have all the info we need to decide if this thing is a spike or not.
        for i = 1:size(spikes, 1)  % for each spike
            if sum(spikes(i,[3 5])) > thresh && sum(spikes(i,[3 5])) > absthresh            % both parts together are bigger than thresh: so have some flexibility in relative sizes
                if spikes(i,2) > spikedur     % spike wave cannot be too sharp: then it is either too small or noise
                    out(end+1,1) = spikes(i,1);         % add timestamp of spike to output list
                else
                    toosharp(end+1) = spikes(i,1);
                end
            else
                toosmall(end+1) = spikes(i,1);
            end
        end

        totalspikes(dd) =  totalspikes(dd) + length(out);  % keep track of total number of spikes so far

        
         %% Re-align spikes to peak of the spikey component
         timeToPeak = [-.1,.15];
         idxToPeak = timeToPeak*srate;
         for i = 1:size(out,1)
            currIdx = out(i,1);
            idxToLook = round(currIdx+idxToPeak(1)):round(currIdx+idxToPeak(2));
            snapshot = HFdata(idxToLook);
            [~,I] = max(abs(snapshot));
            out(i,1) = out(i,1) + idxToPeak(1) + I;
         end
        
        

        if ~isempty(out)
            %error('look\n');
            out(:,2) = dd;
            out(:,1) = out(:,1);
            out(:,1) = out(:,1)+time_points(1)-1;
            allout = [allout; out];
        end
        
        
       
       all_spikes = [all_spikes;allout];
   
    if 1== 0 && dd == 1
        figure
        plot_times = 1:15*srate;
        plot(linspace(0,15,length(plot_times)),data(plot_times))
        hold on
        plot(linspace(0,15,length(plot_times)),HFdata(plot_times))
        error('look\n');
    end
       
     
        
        
    

        
    end
    %fprintf('no');
end


% need to remove duplicates
% need to sort and flip lr

%allout = sortrows(allout);
%allout = fliplr(allout);

gdf = all_spikes;
gdf = unique(gdf,'stable','rows');

if isempty(gdf) == 0
    times = gdf(:,1);
    chs = gdf(:,2);
    [times,I] = sort(times);
    chs = chs(I);
    gdf = [chs,times];
end


%{
% % Plotting
if size(eeg,2) > 5
 offset = linspace(-10000,10000,n_chans);
 for chan =1:n_chans
%     
%     
     plot(eeg(:,chan)+offset(chan),'b'); hold on;
%     
    %{
     idx = find(gdf(:,1)==chan);
     if isempty(idx) == 0
        for i = 1:length(idx)
         tmp  = gdf(idx(i),2);
         text(tmp,eeg(tmp,chan)+offset(chan),'*','color','red','fontsize',35,'fontweight','bold');
        end 
     end
%}
     col = zeros(size(noise,1),3);
     for nn = 1:size(noise,1)
        if noise(nn,chan) == -1, col(nn,:) = [0 0 1]; else, col(nn,:) = [1 0 0]; end
     end
%     
     hold on;
     n_amp = noise(:,chan)*500+offset(chan);
     scatter(linspace(1,size(eeg,1),size(noise,1)),n_amp,100,col,'filled')
%     
 end
% 
 set(gca,'ydir','reverse')
 set(gca,'ylim',[min(offset)-500,max(offset)+500])
 
 fprint('no\n');
end

%}