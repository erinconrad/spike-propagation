function plotNoisyTimes(pt,whichPt,sz)

%% Get paths and load seizure info and channel info
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);

times = pt(whichPt).sz(sz).stats.noise(:,1)+pt(whichPt).sz(sz).runTimes(1,1);
noise = pt(whichPt).sz(sz).stats.noise(:,2);

nanmedian(noise)
noisyIdx = (noise>10*nanmedian(noise));
noisy = noise(noisyIdx);
noisytimes = times(noisyIdx);

szTime = pt(whichPt).sz(sz).onset;

plot(times,noise)
xlabel('Time (s)','FontSize',15);
ylabel('Noise','FontSize',15);
hold on
scatter(noisytimes,noisy,'r')
ylims = get(gca,'ylim');
plot([szTime, szTime],[ylims(1) ylims(2)], 'k');
title(sprintf('Noise over time for %s seizure %d',pt(whichPt).name,sz),'FontSize',15);

outputFile = [pt(whichPt).name,'_sz_',sprintf('%d',sz),'_noise_.png'];

saveas(gcf,[resultsFolder,'plots/',pt(whichPt).name,'/',outputFile])

end