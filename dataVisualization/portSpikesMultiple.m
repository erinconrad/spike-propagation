%{

This function visualizes spike detections for specified times


%}


function portSpikesMultiple(Patient,pt,chIds,tmul,absthresh,startTimes)

%% Parameters
% Try not to get too close to seizure!
redo = 0;


    

% how much time to detect
duration = 600; %default 600, less than this will screw up sensitivity of detector

% how much time to plot
plot_duration = 15; %default 15

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;


dataName = Patient(pt).ieeg_name;
electrodeFile = Patient(pt).electrode_labels;
ptname = Patient(pt).name;

outputFile = [ptname,'_tmul_',sprintf('%d',tmul),'_absthresh_',sprintf('%d',absthresh),'.png'];
outputFolder = [resultsFolder,'spike verification/',Patient(pt).name,'/'];
if exist([outputFolder,outputFile],'file') ~= 0 && redo == 0
    fprintf('Already did %s tmul %d absthresh %d, skipping\n',ptname,tmul,absthresh);
    return
end

for i = 1:length(startTimes)
    time(i).startTime = startTimes(i);
    time(i).times = [startTimes(i),startTimes(i)+duration];
end

whichCh = chIds(1:min(10,length(chIds)));


%% load seizure info and channel info
ptInfo = loadjson(jsonfile);

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;


for i = 1:length(startTimes)
    %% calculate gdf (spike times and locations) and output the data in that time
    fprintf('Detecting spikes\n');
    [time(i).gdf,~,extraoutput] = getSpikeTimes(time(i).times,ptname,dataName,electrodeFile,ptInfo,pwfile,0,0,0,1,0,1,0,tmul,absthresh);
    time(i).values = extraoutput{1};
    time(i).unignoredChLabels = extraoutput{2};
    time(i).plottimes =  [1:size(time(i).values,1)]/fs;
end


%% Plot 
prows = 2;
pcolumns = 2;
colors = {'b','r','g','c','m','b','r','g','c','m'};
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.95, 0.8]);

for s = 1:length(time)
    whichcol = ceil(s/2);
    whichrow = mod(s+1,2)+1;
    subplot('Position',[(whichcol-1)*1/pcolumns (prows-whichrow)*1/prows 1/pcolumns 1/prows])
    range = 0;
    pl = zeros(length(whichCh),1);
    for i = 1:length(whichCh)
        ch = whichCh(i);
        amps = time(s).values(:,ch) - range;
        pl(i) = plot(time(s).plottimes,amps,colors{i});
        hold on

        % find the times of the spikes with the desired channel
        spiketimes = time(s).gdf(time(s).gdf(:,1) == ch,2);

        spikeamp = ones(size(spiketimes,1),1)*max(time(s).values(:,ch))-range;

        scatter(spiketimes,spikeamp,80,colors{i},'filled');
        range = range + max(time(s).values(:,ch)) - min(time(s).values(:,ch));
    end
    legnames = time(s).unignoredChLabels(whichCh);
    legend(pl,legnames,'Location','northeast');
    xlabel('Time (s)');
    xlim([0 plot_duration]);
    set(gca,'YTickLabel',[]);
    text(0.1,0.95,sprintf('%s Time %1.1f, tmul %d, absthresh %d, %d s duration',...
        ptname,time(s).startTime, tmul, absthresh,plot_duration),'units','normalized','FontSize',15);
    set(gca,'fontsize',15);
    xticks = 1:plot_duration-1;
    
    yl = ylim;
    yloc = yl(1)+(yl(2)-yl(1))*0.05;
    for t = 1:length(xticks)
       tt = text(xticks(t),yloc,sprintf('%d s',xticks(t)),'FontSize',15); 
    end


end


saveas(gcf,[outputFolder,outputFile])
close(gcf);

fprintf('\n');



end