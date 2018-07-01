function [gdf,extraOutput] = getSpikesSimple(pt,whichPt,times,whichDetector)

%% Parameters
setChLimits = 1; % Should I toss out spikes that occur across too many channels at the same time
multiChLimit = 0.8; % I will throw out spikes that occur in >80% of channels at the same time
multiChTime = .025;


%% Load file paths, etc.
[~,~,~,~,pwfile] = fileLocations;
dataName = pt(whichPt).ieeg_name;


% Get the indices I want to look at
fs = pt(whichPt).fs;

oldtimes = times;
oldStartEnd = oldtimes*fs;

% Force it to look at least 60 seconds if detector 4
if whichDetector == 4 && times(2)-times(1)<60
    times(2) = times(1)+60;
end

channels = pt(whichPt).channels;
startAndEndIndices = times*fs;
indices = startAndEndIndices(1):startAndEndIndices(2);
tmul = pt(whichPt).tmul;
absthresh = pt(whichPt).absthresh;

%% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices,pwfile);

%% Run spike detector
if whichDetector == 1

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

elseif whichDetector == 3
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

    [gdf] = fspk3(data.values,tmul,absthresh,length(channels),data.fs,window);

    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
        gdf(:,2) = gdf(:,2)/data.fs;
    end

end

%% Toss spikes that occur across too high a percentage of channels at the same time

if setChLimits == 1 && isempty(gdf) == 0
   maxChannels = multiChLimit * length(channels);
   newgdf = tooManyElectrodes(gdf,maxChannels,multiChTime);
   fprintf('Percentage of spikes discarded for being across too many channels: %1.1f\n',...
        (size(gdf,1)-size(newgdf,1))/size(gdf,1)*100)
   gdf = newgdf;
    
end

% limit spikes and values to those in the desired times
values = data.values;
values = values(1:oldStartEnd(2)-oldStartEnd(1),:);
gdf = gdf(gdf(:,2)<oldtimes(2)-oldtimes(1),:);

extraOutput.values = values;


end