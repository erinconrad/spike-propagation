function [sensitivity,accuracy] = spikeChecker(pt,whichPt,chs,...
    spikeTimes,notSpikeTimes,tmul,absthresh,whichDetector,trainOrTest,whichSpikes)

chunkSize = 10;

if whichDetector ==  2
    detect_duration = [-1 599]; %1 seconds before, 599 after
elseif whichDetector == 4
    detect_duration = [-1 59]; %1 seconds before, 59 after 
end
check_duration = 2; % look for 2 seconds after the start of the detection period (which begins 1 s before the spike)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;


dataName = pt(whichPt).ieeg_name;
electrodeFile = pt(whichPt).electrode_labels;
ptname = pt(whichPt).name;


if trainOrTest == 1
    outputFolder = [resultsFolder,'validation/',pt(whichPt).name,'/train/'];
    trainText = 'Training';
elseif trainOrTest == 2
    outputFolder = [resultsFolder,'validation/',pt(whichPt).name,'/test/'];
    trainText = 'Testing';
end
mkdir(outputFolder,sprintf('tmul_%d_absthresh_%d_detector_%d/',tmul,absthresh,whichDetector))
outputFolder = [outputFolder,'/',sprintf('tmul_%d_absthresh_%d_detector_%d/',tmul,absthresh,whichDetector)];


times = [spikeTimes;notSpikeTimes];
allWhichSpikes = [whichSpikes;whichSpikes];

for i = 1:length(times)
    time(i).runTimes = times(i)+detect_duration;
    time(i).times = times(i);
    if i<=length(spikeTimes)
        time(i).isSpike = 1;
    else
        time(i).isSpike = 0;
    end
    time(i).whichSpike = allWhichSpikes(i);
end

%% load seizure info and channel info
ptInfo = loadjson(jsonfile);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% Detect spikes
for i = 1:length(time)
    [time(i).gdf,~,extraoutput] = getSpikeTimes(time(i).runTimes,ptname,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,1,0,tmul,absthresh,whichDetector);
    time(i).values = extraoutput{1};
    
    
    time(i).unignoredChLabels = extraoutput{2};
    time(i).indicesInCheckTime = 1:check_duration*fs;
    time(i).valuesInCheckTime =  time(i).values(time(i).indicesInCheckTime,:);
    
    if isempty(time(i).gdf)==0
        
        time(i).spikesInCheckTime = time(i).gdf(time(i).gdf(:,2)<=check_duration,:);
    else
        time(i).spikesInCheckTime = [];
    end
    time(i).plottimes =  [1:size(time(i).valuesInCheckTime,1)]/fs;  
end

%% Decide if TP, FP, TN, FN
TP = 0;
FP = 0;
TN = 0;
FN = 0;
for i = 1:length(time)
   % if there is a spike there per my gold standard human detection 
   if time(i).isSpike == 1 
      
       
       % If no spikes at all
       if isempty(time(i).spikesInCheckTime) == 1
           time(i).designation = 'FN';
           time(i).color = 'r';
           FN = FN + 1;
           continue
       end
       
       % If I detected at least one spike, in ANY channel
       if isempty(find(ismember(time(i).spikesInCheckTime(:,1),chs),1)) == 0
           time(i).designation = 'TP';
           time(i).color = 'g';
           TP = TP + 1;
       
       % If I didn't detect any spikes at all
       elseif isempty(find(ismember(time(i).spikesInCheckTime(:,1),chs),1)) == 1
           time(i).designation = 'FN';
           time(i).color = 'r';
           FN = FN + 1;
           
       end
       
   % if there is NO spike there per my gold standard human detection
   elseif time(i).isSpike == 0
       
       % If no spikes at all
       if isempty(time(i).spikesInCheckTime) == 1
           time(i).designation = 'TN';
           time(i).color = 'g';
           TN = TN + 1;
           continue
       end
       
       % If I detected at least one spike, in ANY channel
       if isempty(find(ismember(time(i).spikesInCheckTime(:,1),chs),1)) == 0
           time(i).designation = 'FP';
           time(i).color = 'r';
           FP = FP + 1;
       
       % If I didn't detect any spikes at all
       elseif isempty(find(ismember(time(i).spikesInCheckTime(:,1),chs),1)) == 1
           time(i).designation = 'TN';
           time(i).color = 'g';
           TN = TN + 1;
           
       end
       
    % my designation for not-sure
   elseif time(i).isSpike == 2
          time(i).designation = '??';
          time(i).color = 'k';
       
   end
    
end

%% Calculate sensitivity and accuracy

% sensitivity is number of correctly detected spikes divided by number of
% true spikes
sensitivity = TP/(TP+FN);

accuracy = TP/(TP + FP + FN);

%% Do plots
% Divide the times into chunks of 5
nchunks = ceil(length(times)/chunkSize);
for ichunk = 1:nchunks
    spike_idx = (ichunk-1)*chunkSize+1:min(ichunk*chunkSize,length(times));
    outputFile = [ptname,'_tmul_',sprintf('%d',tmul),'_absthresh_',sprintf('%d',absthresh),...
        'chunk_',sprintf('%d',ichunk),'.png'];
    prows = 1;
    pcolumns = length(spike_idx);
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.99, 0.8]);
    for sidx = 1:length(spike_idx)
        s = spike_idx(sidx);
        sp(sidx) = subplot(1,length(spike_idx),sidx);
       
        range = 0;
        pl = zeros(length(chs),1);
        
        if time(s).isSpike == 1
            spikeText = sprintf('Spike %d\n%s\n%s',time(s).whichSpike,time(s).designation,trainText);
        elseif time(s).isSpike == 0
            spikeText = sprintf('Not-a-spike %d\n%s\n%s',time(s).whichSpike,time(s).designation,trainText);
        end

        for i = 1:length(chs)
            ch = chs(i);
            amps = time(s).valuesInCheckTime(:,ch) - range;
            pl(i) = plot(time(s).plottimes,amps,'k');
            hold on

            if isempty(time(s).spikesInCheckTime) == 0
                % find the times of the spikes with the desired channel
                spiketimes = time(s).spikesInCheckTime(time(s).spikesInCheckTime(:,1) == ch,2);

                spikeamp = ones(size(spiketimes,1),1)*max(time(s).valuesInCheckTime(:,ch))-range;

                scatter(spiketimes,spikeamp,80,'k','filled');
            end
            
            range = range + max(time(s).valuesInCheckTime(:,ch)) - min(time(s).valuesInCheckTime(:,ch));
        end

        %legnames = time(s).unignoredChLabels(chs);
        %legend(pl,legnames,'Location','northeast');
        xlabel('Time (s)');
        xlim([0 check_duration]);
        set(gca,'YTickLabel',[]);
      %  text(0.1,0.95,sprintf('%s Time %1.1f, tmul %d, absthresh %d, %d s duration',...
      %      ptname,time(s).startTime, tmul, absthresh,plot_duration),'units','normalized','FontSize',15);

        %xticks = 1:check_duration-1;

        %yl = ylim;
        %yloc = yl(1)+(yl(2)-yl(1))*0.05;
        %for t = 1:length(xticks)
        %   tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
        %end

        %fprintf('check\n');
        text(0.5,0.95,spikeText,'units','normalized','FontSize',20,...
            'color',time(s).color,'HorizontalAlignment','center');
        set(gca,'fontsize',15);

    end
    
    for sidx = 1:length(spike_idx)
       
        whichcol = sidx;
        whichrow = 1;
        
        set(sp(sidx),'Position',[(whichcol-1)*1/pcolumns (prows-whichrow)*1/prows 1/pcolumns 1/prows])
    end

    saveas(gcf,[outputFolder,outputFile])

end

%fprintf('test\n');

end