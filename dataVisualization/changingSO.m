function changingSO(P,pt,sz,whichToPlot)    

%{ 
This function plots sequence frequency, spatial organization, and the ratio
of the two over time surrounding a seizure


%}

scale = 3600;
allSO = [];
allSF = [];
times = [];
for i = 1:length(P(pt).sz(sz).tblock)
    allSO = [allSO,P(pt).sz(sz).tblock(i).spatialOrg];
    times = [times,(P(pt).sz(sz).tblock(i).times(1)+...
        P(pt).sz(sz).tblock(i).times(2))/2+P(pt).sz(sz).runTimes(1)];
    
    allSF = [allSF,length(P(pt).sz(sz).tblock(i).sIdx)];
end

allRat =  allSO./allSF;

szTimes = P(pt).sz(sz).szTimes;
insz = zeros(size(times));

for i = 1:length(P(pt).sz(sz).tblock)
   bothtimes = P(pt).sz(sz).tblock(i).times+P(pt).sz(sz).runTimes(1);
   if szTimes(2) > bothtimes(1) && szTimes(1) < bothtimes(2)
       insz(i) = 1;
   end
    
end

tsz = times(insz==1);

if whichToPlot == 1
    toPlot = allSF;
    ptitle = 'Sequence frequency';
    pname = 'SF';
elseif whichToPlot == 2
    toPlot = allSO;
    ptitle = 'Spatial organization';
    pname = 'SO';
elseif whichToPlot == 3
    toPlot = allRat;
    ptitle = 'Ratio of spatial organization to spike frequency';
    pname = 'rat';
end

plot(times/scale,toPlot,'color','k','linewidth',2);
hold on
ylims = get(gca,'ylim');
patch([min(tsz) max(tsz) max(tsz) min(tsz)]/scale,[ylims(1) ylims(1) ylims(2) ylims(2)],...
    [0.6 0.6 0.6],'Edgecolor','none');
plot(times/scale,toPlot,'color','k','linewidth',2);
plot([szTimes(1) szTimes(1)]/scale,get(gca,'ylim'),'--','color','k','linewidth',3);
ylabel(sprintf('%s',ptitle));
xlabel('Time (h)');
title([sprintf('%s surrounding seizure for patient HUP%d',ptitle,pt)]);
set(gca,'FontSize',15)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
outputFile = ['HUP_',sprintf('%d',pt),'_','_sz_',sprintf('%d',sz)...
    ,sprintf('%s',pname),'.png'];
saveas(gcf,[resultsFolder,outputFile])

end