%{
This is my alternate function for detecting EKG spikes. This is basically
the usual detector without looking for an aftergoing slow wave

%}

function spikes = detectEKGSpikes(eeg,fs)

%% Parameters
window = 60 * fs; % window over which to assess baseline
tmul = 10; % threshold multiplier
spkdur = 220;
fr     = 20; 
lfr    = 7;  

%% Initialize stuff
n_chans = size(eeg,2);
spikes = [];
num_segs = ceil(size(eeg,1)/window);

%% Loop through EKG channels and times and detect spikes
for dd = 1:n_chans
   for tt = 1:num_segs 
       time_points(1) = (tt-1)*window+1;
       time_points(2) = min(window*tt,size(eeg,1));
       if tt == num_segs-1
           time_points(2) = size(eeg,1);
       end
       
       % Get the data chunk for that one minute and that channel
       data = eeg(time_points(1):time_points(2),dd);
       
       % Get the threshold that the spike needs to exceed
        lthresh = mean(abs(data));  % this is the smallest the initial part of the spike can be
        thresh  = lthresh*tmul;     % this is the final threshold we want to impose
        sthresh = lthresh*tmul/3;   % this is the first run threshold
      
       % perform high pass filter
       fndata   = eegfilt(data, 1, 'hp',fs);
       HFdata    = eegfilt(fndata, fr, 'lp',fs);
       
       % Find peaks
       [spp,spv] = FindPeaks(HFdata);
       
       % check that they're close enough together
       idx      = find(diff(spp) <= spkdur); 
       startdx  = spp(idx);
       startdx1 = spp(idx+1);
       
       % check the amplitude of the waves of appropriate duration
        for i = 1:length(startdx)
            spkmintic = spv(find(spv > startdx(i) & spv < startdx1(i)));  % find the valley that is between the two peaks

            if HFdata(startdx1(i)) - HFdata(spkmintic) > sthresh & HFdata(startdx(i)) - HFdata(spkmintic) > lthresh 
                if HFdata(startdx1(i)) - HFdata(spkmintic) > thresh
                 spikes(end+1) = spkmintic+time_points(1); 
                end % add timestamp to the spike list
            end

        end
    
    
   end 
end





end