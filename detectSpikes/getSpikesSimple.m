function [gdf,extraOutput] = getSpikesSimple(pt,whichPt,times,whichDetector,thresh)

% This is my basic spike detector file which can call one of a number of
% specific detectors

% Need to uncomment lines 192 and 200 to get vanleer back

%% Parameters
setChLimits = 1; % Should I toss out spikes that occur across too many channels at the same time
multiChLimit = 0.8; % I will throw out spikes that occur in >80% of channels at the same time
multiChTime = .025;
vtime = [-0.002,0.048]; % Do I use this?

%% Load file paths, etc.
[~,~,~,~,pwfile] = fileLocations;
dataName = pt(whichPt).ieeg_name;

% If no thresh, use a default
if exist('thresh','var') == 0
    error('Warning, no thresholds entered, using tmul 13 and absthresh 300\n');
    thresh.tmul = 13;
    thresh.absthresh = 300;
end

% Get the indices I want to look at
fs = pt(whichPt).fs;
oldtimes = times;
oldStartEnd = oldtimes*fs;

% chLocs
chLocs = pt(whichPt).electrodeData.locs;

% Force it to look at least 60 seconds if detector 4
if whichDetector == 4 && times(2)-times(1)<60
    fprintf('Warning, need at least 60 seconds to use detector 4. Running with 60 s of data.');
    times(2) = times(1)+60;
end

channels = pt(whichPt).channels;
startAndEndIndices = times*fs;
indices = startAndEndIndices(1):startAndEndIndices(2);
tmul = thresh.tmul;
absthresh = thresh.absthresh;

%% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices,pwfile);

% remove nans
data.values(isnan(data.values)) = 0;

%% Notch filter to remove 60 Hz noise
f = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',fs);
           
for i = 1:size(data.values,2)
   data.values(:,i) = filtfilt(f,data.values(:,i));   
end


%% Run spike detector
if whichDetector == 1
    
    fprintf('Warning, why are you using Detector 1?\n');

    % This is the spike detector from Janca et al 2014, edited by me as
    % above
    [out,~,~,~,~,~] = spike_detector_Erin(data.values,data.fs);
    % reorder spikes by time
    [timeSort,I] = sort(out.pos);
    chanSort = out.chan(I);
    % make gdf
    if isempty(out.pos) == 1
        fprintf('Warning: No spikes detected\n');
    else
        %fprintf('Detected %d spikes\n',length(out.pos));
        gdf = [chanSort,timeSort];
    end

elseif whichDetector == 2
    
    fprintf('Warning, why are you using Detector 2?\n');
     addpath('./SamCode');
    % This calls the Bermudez detector
    % 
    % I have not edited this at all at this point.
    gdf = fspk2(data.values,tmul,absthresh,length(channels),data.fs);
    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
         gdf(:,2) = gdf(:,2)/data.fs;

        out.pos = gdf(:,2); out.chan = gdf(:,1);
    end
    
    noise = [];
    removed = [];

elseif whichDetector == 3
    fprintf('Warning, why are you using Detector 3?\n');
    
    %this is the unedited Janca detector
    [out,~,~,~,~,~] = spike_detector_hilbert_v16_nodownsample(data.values,data.fs,'-h 60');
    % reorder spikes by time
    [timeSort,I] = sort(out.pos);
    chanSort = out.chan(I);

    % make gdf
    if isempty(out.pos) == 1
        fprintf('Warning: No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',length(out.pos));
        gdf = [chanSort,timeSort];
    end

elseif whichDetector == 4
    % this is my edited version of Sam's spike detector, which instead of
    % measuring the relative threshold based on the absolute value of the
    % entire data, just looks at a minute surrounding the potential spike
    window = 60*data.fs;

    [gdf,noise,removed] = fspk3(data.values,tmul,absthresh,...
        length(channels),data.fs,window);

    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
        gdf(:,2) = gdf(:,2)/data.fs;
    end
    
    if isempty(removed) == 0
        removed(:,2:3) = removed(:,2:3)/data.fs;
    end
    
    
elseif whichDetector == 5
    % just like fspk3 (detector 4) except I am changing fr and I have the ability to
    % change parameters based on if depth
    window = 60*data.fs;

    [gdf,noise,removed] = fspk4(data.values,tmul,absthresh,...
        length(channels),data.fs,window,pt(whichPt).electrodeData.electrodes);

    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
        gdf(:,2) = gdf(:,2)/data.fs;
    end
    
    if isempty(removed) == 0
        removed(:,2:3) = removed(:,2:3)/data.fs;
    end

end


%% Re-align the spike to be the peak
values = data.values;
timeToPeak = [-.1,.15]; % Where to look for the peak
idxToPeak = timeToPeak*fs;

new_gdf = gdf;

% Loop through spikes
for i = 1:size(gdf,1)
    snapshot = values(max(gdf(i,2)*fs+idxToPeak(1),1):gdf(i,2)*fs+idxToPeak(2),gdf(i,1));
    [~,I] = max(snapshot);
    new_gdf(i,2) = gdf(i,2) + timeToPeak(1) + I/fs;
end

gdf = new_gdf;

% Re-sort spike times
if isempty(gdf) == 0
    times_t = gdf(:,2);
    chs = gdf(:,1);
    [times_t,I] = sort(times_t);
    chs = chs(I);
    gdf = [chs,times_t];
end



% test plot for re-aligned spikes
%{
for whichSp = 10:min(size(gdf,1),30)
toplot = values((gdf(whichSp,2)*fs-2*fs:gdf(whichSp,2)*fs+13*fs),gdf(whichSp,1));
plot(linspace(-2,13,length(toplot)),toplot);
hold on
scatter(0,values(round(gdf(whichSp,2)*fs),gdf(whichSp,1)));
scatter(new_gdf(whichSp,2)-gdf(whichSp,2),values(round(new_gdf(whichSp,2)*fs),new_gdf(whichSp,1)))
hold off
pause
end
%}


%% Toss spikes that occur across too high a percentage of channels at the same time
if setChLimits == 1 && isempty(gdf) == 0
   maxChannels = multiChLimit * length(channels);
   newgdf = tooManyElectrodes(gdf,maxChannels,multiChTime);
   fprintf('Percentage of spikes discarded for being across too many channels: %1.1f\n',...
        (size(gdf,1)-size(newgdf,1))/size(gdf,1)*100)
   gdf = newgdf;
    
end

%% limit spikes and values to those in the desired times
% This is useful if I had to force lengthen the search time for detector 4
% if I asked to look at less than 60 seconds of data
values = values(1:oldStartEnd(2)-oldStartEnd(1),:);
if isempty(gdf) == 0
    gdf = gdf(gdf(:,2)<oldtimes(2)-oldtimes(1),:);
end


extraOutput.values = values;

if isempty(gdf) == 1
    vanleer = [];
else
    vanleer = [];% vMakeSegments(gdf,data.values,fs,vtime,chLocs);
end



% Re-align gdf times to be the actual times
if isempty(gdf) == 0
    gdf(:,2) = gdf(:,2) + times(1);
    vanleer = [];
    %vanleer.spikeTimes = vanleer.spikeTimes + times(1);
end

if isempty(removed) == 0
    removed(:,2:3) = removed(:,2:3) + times(1);
end


extraOutput.vanleer = vanleer;
extraOutput.removed = removed;

[~,noisychs] = find(noise == 1);
extraOutput.noise.noisychs = pt(whichPt).electrodeData.unignoredChs(noisychs);
extraOutput.noise.noisematrix = noise;


end