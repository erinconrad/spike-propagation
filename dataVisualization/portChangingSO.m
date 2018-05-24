function portChangingSO(P,pt,sz,whichToPlot)    

%{ 
This function plots sequence frequency, spatial organization, and the ratio
of the two over time surrounding a seizure


%}

scale = 3600;
allSO = [];
allSF = [];
times = [];
allCI = [];
for i = 1:length(P(pt).sz(sz).blockRL)
    allSO = [allSO,P(pt).sz(sz).blockRL(i).spatialOrg];
    times = [times,(P(pt).sz(sz).blockRL(i).times(1)+...
        P(pt).sz(sz).blockRL(i).times(2))/2+P(pt).sz(sz).runTimes(1)];
    
    allSF = [allSF,length(P(pt).sz(sz).blockRL(i).sIdx)];%/...
  %      (P(pt).sz(sz).blockRL(i).times(2)-P(pt).sz(sz).blockRL(i).times(1))];
  
    allCI = [allCI;P(pt).sz(sz).blockRL(i).CI95];
end

allRat =  allSO./allSF;

szTimes = [P(pt).sz(sz).onset,P(pt).sz(sz).offset];
insz = zeros(size(times));

for i = 1:length(P(pt).sz(sz).blockRL)
   bothtimes = P(pt).sz(sz).blockRL(i).times+P(pt).sz(sz).runTimes(1);
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

if whichToPlot == 2
   plot(times/scale,allCI(:,1),'--')
   plot(times/scale,allCI(:,2),'--')
end

ylabel(sprintf('%s',ptitle));
xlabel('Time (h)');
title([sprintf('%s',ptitle),' surrounding seizure for ',sprintf('%s',P(pt).name)]);
set(gca,'FontSize',15)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
outputFile = [sprintf('%s',P(pt).name),'_','_sz_',sprintf('%d',sz)...
    ,sprintf('%s',pname),'.png'];
saveas(gcf,[resultsFolder,'plots/',outputFile])

end