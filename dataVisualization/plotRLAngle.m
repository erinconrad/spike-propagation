function plotRLAngle(P,pt,sz)
%{ 

This script plots recruitment latencies over the electrodes



%}

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
ptname = P(pt).name;
dotsize = 100;
biggerdotsize = 200;


% Get the channel locations for the patient
chLocs = P(pt).sz(sz).data.xyChan(:,2:4);

earlychans_interic = P(pt).sz(sz).ic.earlychans_interic;
earlychans_ic = P(pt).sz(sz).ic.earlychans_ic;
latechans_interic = P(pt).sz(sz).ic.latechans_interic;
latechans_ic = P(pt).sz(sz).ic.latechans_ic;
direction_interic = P(pt).sz(sz).ic.direction_interic;
direction_ic = P(pt).sz(sz).ic.direction_ic;



% Show all electrode locations, making them unfilled and dark
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.4, 0.5]);
rLat = P(pt).sz(sz).ic.rl_interictal;
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,rLat,'filled','MarkerEdgeColor','k')
hold on
plot3([earlychans_interic(1),latechans_interic(1)],[earlychans_interic(2),latechans_interic(2)],...
    [earlychans_interic(3),latechans_interic(3)],'k','LineWidth',3);
e=scatter3(earlychans_interic(1),earlychans_interic(2),earlychans_interic(3),...
    biggerdotsize,'gs','filled','MarkerEdgeColor','k','LineWidth',1);
l=scatter3(latechans_interic(1),latechans_interic(2),latechans_interic(3),...
    biggerdotsize,'rd','filled','MarkerEdgeColor','k','LineWidth',1);

title('Electrode activation delay for interictal sequences');
legend([e,l],{'Mean position early electrodes','Mean position late electrodes'},'location','northwest',...
    'Fontsize',15);
c = colorbar;
%colormap jet
%set(gca,'Position',[0.1300 0.1100 0.80 0.8]);
set(c,'FontSize',15);
c.Label.String = 'seconds';
%c.Label.FontSize = 18;
grid off
axis off
set(gca,'FontSize',15)

outputFile = [ptname,'_','_sz_',sprintf('%d',sz),...
    '_interictal_recruitmentLatency','.eps'];

%saveas(gcf,[resultsFolder,'plots/',ptname,'/',outputFile])
print([resultsFolder,'plots/',ptname,'/',outputFile],'-depsc');


figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.4, 0.4, 0.5]);
rLat = P(pt).sz(sz).ic.rl_ictal;
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,rLat,'filled','MarkerEdgeColor','k')
hold on
plot3([earlychans_ic(1),latechans_ic(1)],[earlychans_ic(2),latechans_ic(2)],...
    [earlychans_ic(3),latechans_ic(3)],'k','LineWidth',3);
e=scatter3(earlychans_ic(1),earlychans_ic(2),earlychans_ic(3),...
    biggerdotsize,'gs','filled','MarkerEdgeColor','k','LineWidth',1);
l=scatter3(latechans_ic(1),latechans_ic(2),latechans_ic(3),...
    biggerdotsize,'rd','filled','MarkerEdgeColor','k','LineWidth',1);

title('Electrode activation delay for ictal sequences');
legend([e,l],{'Mean position early electrodes','Mean position late electrodes'},'location','northwest',...
    'Fontsize',15);
c=colorbar;
set(c,'FontSize',15);
c.Label.String = 'seconds';
%c.Label.FontSize = 18;

%colormap jet
grid off
axis off
set(gca,'FontSize',15)


outputFile = [ptname,'_','_sz_',sprintf('%d',sz),...
    '_ictal_recruitmentLatency','.eps'];

%saveas(gcf,[resultsFolder,'plots/',ptname,'/',outputFile])
print([resultsFolder,'plots/',ptname,'/',outputFile],'-depsc');



end