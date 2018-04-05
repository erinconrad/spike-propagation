function plotRecruitmentLatency(P,sz,block)
%{ 

This script plots recruitment latencies over the electrodes

%}
ptname = 'HUP80';
dotsize = 70;
avgRecruitmentLat = P(80).sz(sz).block(block).data.avgRecruitmentLat;

% Get the channel locations for the patient
chLocs = P(80).sz(sz).block(block).data.xyChan(:,2:4);
colormin = min(avgRecruitmentLat);
colormax = max(avgRecruitmentLat);

% Show all electrode locations, making them unfilled and dark
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,'k')

hold on

for i = 1:size(chLocs,1)
    rLat= avgRecruitmentLat(i);
    if isnan(rLat) == 0
        scatter3(chLocs(i,1),chLocs(i,2),chLocs(i,3),dotsize,rLat,'filled')
        
    end
    
end

colormap jet
colorbar


outputFile = [ptname,'_','_sz_',sprintf('%d',sz),'_block_',...
    sprintf('%d',block),'_recruitmentLatency','.png'];

saveas(gcf,[resultsFolder,outputFile])


end