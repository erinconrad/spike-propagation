function plotRecruitmentLatency(P,pt,sz,block)
%{ 

This script plots recruitment latencies over the electrodes



%}

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
ptname = P(pt).name;
dotsize = 100;
%avgRecruitmentLat = P(pt).sz(sz).blockRL(block).avgRecruitmentLat;
avgRecruitmentLat = P(pt).sz(sz).ic.rl_interictal;
%avgRecruitmentLat = P(pt).sz(sz).MI.rl;

% Get the channel locations for the patient
chLocs = P(pt).sz(sz).data.xyChan(:,2:4);
colormin = min(avgRecruitmentLat);
colormax = max(avgRecruitmentLat);

% Show all electrode locations, making them unfilled and dark
figure
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,'k')

hold on

for i = 1:size(chLocs,1)
    rLat= avgRecruitmentLat(i);
    if isnan(rLat) == 0
        scatter3(chLocs(i,1),chLocs(i,2),chLocs(i,3),dotsize,rLat,'filled')
        
    end
    
end
title(sprintf('Average recruitment latency across channels for %s block %d',ptname, block));
colormap jet
colorbar
grid off
axis off
set(gca,'FontSize',15)

outputFile = [ptname,'_','_sz_',sprintf('%d',sz),'_block_',...
    sprintf('%d',block),'_recruitmentLatency','.png'];

saveas(gcf,[resultsFolder,'plots/',ptname,'/',outputFile])


end