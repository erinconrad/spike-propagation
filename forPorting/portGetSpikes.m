function [gdf,vanleer,noise] = portGetSpikes(desiredTimes,dataName,...
    channels,pwfile,tmul,absthresh,whichDetector,fs)

%{
This is my primary function to detect spikes and output them to a gdf 
with spike times and locations. gdf is an nx2 array, where n is the number
of spikes, the first column has the channel location of the spike and the
2nd column has the time (in seconds) of the spike
%}

%% Parameters
% whichDetector %1 = modified Janca detector, 2 = Bermudez detector, 3 = orig Janca, 4 = modified Bermudez detector
setChLimits = 1; % Should I toss out spikes that occur across too many channels at the same time
multiChLimit = 0.8; % I will throw out spikes that occur in >80% of channels at the same time
multiChTime = .025; % The time period over which spikes need to occur across multiple channels to toss

% segment time (50 ms segments going from 2 ms before the spike detection
% to 48 ms after the spike detection)
vtime = [-0.002,0.048];


%% Prep what data I want to look at

% Get the indices I want to look at
startAndEndIndices = desiredTimes*fs;
indices = startAndEndIndices(1):startAndEndIndices(2);

% get the data from those indices and channels (ignoring ignored channels)
data = getiEEGData(dataName,channels,indices,pwfile);

noise = [];

%% Get bad times
%{
bad = getBadTimes(data);
bad.empty = bad.empty/data.fs;
bad.noise(:,1) = bad.noise(:,1)/data.fs;
%}

%% Run spike detector
if whichDetector == 1

    % This is the spike detector from Janca et al 2014, edited by me as
    % above
    [out,~,~,~,~,~] = spike_detector_Erin(data.values,fs);

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
    
    % I have not edited this at all at this point.
    gdf = fspk2(data.values,tmul,absthresh,length(channels),fs);


    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
        gdf(:,2) = gdf(:,2)/data.fs;
    end
   

elseif whichDetector == 3
    %this is the unedited Janca detector

    [out,~,~,~,~,~] = spike_detector_hilbert_v16_nodownsample(data.values,fs,'-h 60');

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
    window = 60*fs;
    
    [gdf,noise] = fspk3(data.values,tmul,absthresh,length(channels),fs,window);
    
    %{
    [r,~] = (find(noise==1));
    r=unique(r);
    n_times = r*60+desiredTimes(1);
    fprintf('%d is a noisy time \n',n_times);
    %}
    
    if isempty(gdf) == 1
        fprintf('No spikes detected\n');
    else
        fprintf('Detected %d spikes\n',size(gdf,1));
         % put it in seconds
        gdf(:,2) = gdf(:,2)/fs;
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

%{
if setChLimits == 1 && isempty(gdf) == 0
    newgdf =  gdf;
    tooManyChs = [];
    i = 1; % start with i = 1 (the first spike)
    % Loop through the spikes (we are going to variably move through the spikes)
    while 1
        scount = 1;

        % if i is the last spike, break
        if i == size(newgdf,1)
            break
        end

       % for each spike, loop through subsequent spikes to see how many
       % there are across multiple channels occuring at the same time
       for j = i+1:size(newgdf,1)

           % Get the time difference between the first spike and the last
           % spike
           tdiff = newgdf(j,2)-newgdf(i,2);

           % If the difference is small enough, increase the count
           if tdiff < multiChTime
               scount = scount + 1;
           end

           % if the total number of channels spiking in this very close
           % proximity is >80% of the total number of channels
           if scount > multiChLimit*nchan

                % Remove these spikes
                tooManyChs = [tooManyChs;newgdf(i:j,:)];
                newgdf(i:j,:) = [];

                % Break the inner for loop and move to the next spike
                i = j;
                break
           end

       end

       % advance the index of the spike
       i = i+1;
    end
    fprintf('%d total spikes detected\n',size(gdf,1));
    fprintf('Percentage of spikes discarded for being across too many channels: %1.1f\n',...
        size(tooManyChs,1)/size(gdf,1)*100)
    fprintf('%d spikes remain\n',size(newgdf,1));
    gdf = newgdf;
end

%}


if isempty(gdf) == 1
    vanleer = [];
else
    vanleer = vMakeSegments(gdf,data.values,fs,vtime);
end


end