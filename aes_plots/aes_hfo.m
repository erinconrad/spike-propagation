function aes_hfo(pt)

whichPt = 3;
s = 140;
surroundtime = 5;

% Save file location
[~,~,scriptFolder,resultsFolder,pwfile,other] = fileLocations;
destFolder = [resultsFolder,'clustering/plots/'];
mkdir(destFolder);
addpath(genpath(scriptFolder))
addpath(genpath(other.ieeg))

seq = pt(whichPt).seq_matrix(:,s);

% Get the non nan times for the sequence
spike_chs = find(isnan(seq) == 0);
spike_times = seq(find(isnan(seq) == 0));

% sort them from earliest to latest times in the sequence
unsorted_times = spike_times;
[spike_times,I] = sort(spike_times);

unsorted_chs = spike_chs;
spike_chs = spike_chs(I);

times = [spike_times(1)-surroundtime,spike_times(1)+surroundtime];

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
dataName = pt(whichPt).ieeg_name;
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% Get the data for these times
thresh.tmul = pt(whichPt).thresh.tmul;
thresh.absthresh = pt(whichPt).thresh.absthresh;
[gdf,extraoutput] = getSpikesSimple(pt,whichPt,times,6,thresh,0);
values = extraoutput.values;
plottimes =  [1:size(values,1)]/fs;

%% Remove first second
values = values(fs:end-fs,:);
plottimes = plottimes(fs:end-fs);

%% Plot

figure
toAdd = 0;
y_locs = zeros(length(spike_chs),1);
maxy = 0; miny = 0;
for i = 1:length(spike_chs)
    ch = unsorted_chs(i);
    amps = values(:,ch) + toAdd;
    pl(i) = plot(plottimes,amps,'color','k','LineWidth',2);
    hold on
    if max(amps) > maxy, maxy = max(amps); end
    if min(amps) < miny, miny = min(amps); end
    s_time = unsorted_times(i);
    y_locs(i) = toAdd;
    if i ~= length(spike_chs)
       toAdd = toAdd + min(values(:,ch)) - max(values(:,spike_chs(i+1)))-500; 
    end
end
yticklabels([])
xticks([plottimes(1) plottimes(end)])
xticklabels({'0 s',sprintf('%d s',surroundtime*2-2)})
ylim([miny-100 maxy+100])
set(gca,'FontSize',30)

end