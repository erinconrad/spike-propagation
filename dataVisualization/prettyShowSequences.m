function prettyShowSequences(pt,whichPt,sz,s)

doColor = 0;
surroundtime = 1;
dataName = pt(whichPt).ieeg_name;
ptname = pt(whichPt).name;

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);

destFolder = [resultsFolder,'pretty_plots/Fig1/'];


if isempty(sz) ==1
    seq = pt(whichPt).seq_matrix(:,s);
else
    seq = pt(whichPt).sz(sz).seq_matrix(:,s);
end


% Get the non nan times for the sequence
spike_chs = find(isnan(seq) == 0);
spike_times = seq(find(isnan(seq) == 0));

% get colors
if doColor==1
    cmp = jet(length(spike_times));
else
    cmp = zeros(length(spike_times),3);
end

% sort them from earliest to latest times in the sequence
unsorted_times = spike_times;
[spike_times,I] = sort(spike_times);

unsorted_chs = spike_chs;
spike_chs = spike_chs(I);

times = [spike_times(1)-surroundtime,spike_times(1)+surroundtime];

%% Load EEG data info
% calling this with 0 and 0 means I will just get basic info like sampling
% rate and channel labels
data = getiEEGData(dataName,0,0,pwfile);  
fs = data.fs;

%% Get the data for these times
thresh.tmul = pt(whichPt).thresh.tmul;
thresh.absthresh = pt(whichPt).thresh.absthresh;
[gdf,extraoutput] = getSpikesSimple(pt,whichPt,times,6,thresh,0);
values = extraoutput.values;
unignoredChLabels = pt(whichPt).electrodeData.unignoredChs;
plottimes =  [1:size(values,1)]/fs;

sp_vals = zeros(length(spike_times),1);
for i = 1:length(spike_times)
    sp_vals(i) = values((spike_times(i)-times(1))*fs,spike_chs(i));
end

sp_vals_us = zeros(length(spike_times),1);
for i = 1:length(spike_times)
    sp_vals_us(i) = values((unsorted_times(i)-times(1))*fs,unsorted_chs(i));
end

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
    if isequal(cmp(i,:),[1 1 0]) == 1
        cmp(i,:) = [0.8 0.8 0];
    end
    pl(i) = plot(plottimes,amps,'color',cmp(i,:),'LineWidth',2);
    hold on
    if max(amps) > maxy, maxy = max(amps); end
    if min(amps) < miny, miny = min(amps); end
    s_time = unsorted_times(i);
    scatter(s_time-times(1),sp_vals_us(i)+toAdd,80,'k','filled');
    y_locs(i) = toAdd;
    if i ~= length(spike_chs)
       toAdd = toAdd + min(values(:,ch)) - max(values(:,spike_chs(i+1)))-500; 
    end
end
yticks(sort(y_locs))
yticklabels(arrayfun(@(n) sprintf('%d',n),length(spike_chs):-1:1,'UniformOutput',false))
xticks([plottimes(1) plottimes(end)])
xticklabels({'0 s',sprintf('%d s',surroundtime*2-2)})
ylim([miny-100 maxy+100])
set(gca,'FontSize',20)
if doColor==1
    print(gcf,[destFolder,'unsortedColor',sprintf('%d',surroundtime*2-2)],'-depsc');
    eps2pdf([destFolder,'unsortedColor',sprintf('%d',surroundtime*2-2),'.eps'])
else
    print(gcf,[destFolder,'unsortedBW',sprintf('%d',surroundtime*2-2)],'-depsc');
    eps2pdf([destFolder,'unsortedBW',sprintf('%d',surroundtime*2-2),'.eps'])
end
close(gcf)

figure
toAdd = 0;
y_locs = zeros(length(spike_chs),1);
maxy = 0; miny = 0;
for i = 1:length(spike_chs)
    ch = spike_chs(i);
    amps = values(:,ch) + toAdd;
    if isequal(cmp(i,:),[1 1 0]) == 1
        cmp(i,:) = [0.8 0.8 0];
    end
    pl(i) = plot(plottimes,amps,'color',cmp(I(i),:),'LineWidth',2);
    hold on
    if max(amps) > maxy, maxy = max(amps); end
    if min(amps) < miny, miny = min(amps); end
    s_time = spike_times(i);
    scatter(s_time-times(1),sp_vals(i)+toAdd,80,'k','filled');
    y_locs(i) = toAdd;
    if i ~= length(spike_chs)
       toAdd = toAdd + min(values(:,ch)) - max(values(:,spike_chs(i+1)))-500; 
    end
end
yticks(sort(y_locs))
yticklabels(arrayfun(@(n) sprintf('%d',n),flipud(I),'UniformOutput',false))
xticks([plottimes(1) plottimes(end)])
xticklabels({'0 s',sprintf('%d s',surroundtime*2-2)})
set(gca,'FontSize',20)
ylim([miny-100 maxy+100])
if doColor == 1
    print(gcf,[destFolder,'sortedColor',sprintf('%d',surroundtime*2-2)],'-depsc');
    eps2pdf([destFolder,'sortedColor',sprintf('%d',surroundtime*2-2),'.eps'])
else
    print(gcf,[destFolder,'sortedBW',sprintf('%d',surroundtime*2-2)],'-depsc');
    eps2pdf([destFolder,'sortedBW',sprintf('%d',surroundtime*2-2),'.eps'])
end
close(gcf)

end