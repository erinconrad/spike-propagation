%% fspk2 By: Camilo Bermudez 7/31/13

function [gdf] = fspk2(eeg,tmul,absthresh,n_chans,srate)
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
    tmul=13;
end

if ~exist('absthresh')
    absthresh=300;
end

% Initialize parameters
rate   = srate;
chan   = 1:n_chans;
fr     = 20;                 % high pass freq: probably don't have to change
lfr    = 7;                  % low pass freq: probably don't have to change
spkdur = 220;                % spike duration must be less than this
spkdur = spkdur*rate/1000;   % convert to points;

%sharpdur = 200;  % duration of sharp wave
%sharpdur = sharpdur*rate/1000; % convert to points;

aftdur = 150;                % afterhyperpolarization wave must be longer than this
aftdur = aftdur*rate/1000;   % convert to points;

% Not clear on this section...
allout      = [];
overlap     = rate;
totalspikes = zeros(1, length(chan));
inc         = 0;

% Read in eeg data
alldata     = eeg;            % timepnts x chans

% Iterate channels and detect spikes
for dd = 1:n_chans
    
    out     = [];
    data    = alldata(:,dd);
    
    lthresh = mean(abs(data));  % this is the smallest the initial part of the spike can be
    thresh  = lthresh*tmul;     % this is the final threshold we want to impose
    sthresh = lthresh*tmul/3;   % this is the first run threshold
    
    spikes   = [];
    shspikes = [];
    fndata   = eegfilt(data, 1, 'hp',srate);
    
    % first look at the high frequency data for the 'spike' component
    HFdata    = eegfilt(fndata, fr, 'lp',srate);
    [spp,spv] = FindPeaks(HFdata);
    
    idx      = find(diff(spp) <= spkdur);       % find the durations less than or equal to that of a spike
    startdx  = spp(idx);
    startdx1 = spp(idx+1);
    
    % check the amplitude of the waves of appropriate duration
    for i = 1:length(startdx)
        spkmintic = spv(find(spv > startdx(i) & spv < startdx1(i)));  % find the valley that is between the two peaks
        
        if HFdata(startdx1(i)) - HFdata(spkmintic) > sthresh & HFdata(startdx(i)) - HFdata(spkmintic) > lthresh  %#ok<AND2> % see if the peaks are big enough
            spikes(end+1,1) = spkmintic;                                  % add timestamp to the spike list
            spikes(end,2)   = (startdx1(i)-startdx(i))*1000/rate;         % add spike duration to list
            spikes(end,3)   = HFdata(startdx1(i)) - HFdata(spkmintic);    % add spike amplitude to list
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
            else                                         % might be a spike so get the amplitude of the slow wave
                spikes(i,4) = (a(2)-a(1))*1000/rate;       % add duration of afhp to the list
                b = hyperv(find(hyperv > a(1) & hyperv < a(2))); % this is the valley
                spikes(i,5) = LFdata(a(1)) - LFdata(b);  % this is the amplitude of the afhp
                if a(1) == olda                         % if this has the same afterhyperpolarization peak as the prev
                    dellist(end+1) = i-1;               % spike then the prev spike should be deleted
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
    for i = 1:size(spikes, 1);  % for each spike
        if sum(spikes(i,[3 5])) > thresh && sum(spikes(i,[3 5])) > absthresh            % both parts together are bigger than thresh: so have some flexibility in relative sizes
            if spikes(i,2) > 20                     % spike wave cannot be too sharp: then it is either too small or noise
                out(end+1,1) = spikes(i,1);         % add timestamp of spike to output list
            else
                toosharp(end+1) = spikes(i,1);
            end
        else
            toosmall(end+1) = spikes(i,1);
        end
    end
    
    totalspikes(dd) =  totalspikes(dd) + length(out);  % keep track of total number of spikes so far
    
    
    if ~isempty(out)
        out(:,2) = dd;
        out(:,1) = out(:,1);
        allout = [allout; out];
    end
end


% need to remove duplicates
% need to sort and flip lr

allout = sortrows(allout);
allout = fliplr(allout);

gdf = allout;


% % Plotting
% offset = linspace(-10000,10000,n_chans);
% for chan = 1:n_chans
%     
%     
%     plot(eeg(chan,:)+offset(chan),'b'); hold on;
%     
%     idx = find(gdf(:,1)==chan);
%     for i = idx
%         tmp  = gdf(i,2);
%         text(tmp,eeg(chan,tmp)+offset(chan),'*','color','red','fontsize',35,'fontweight','bold');
%     end
%     
%     hold on;
%     
% end
% 
% set(gca,'ydir','reverse')
% set(gca,'ylim',[min(offset)-500,max(offset)+500])

