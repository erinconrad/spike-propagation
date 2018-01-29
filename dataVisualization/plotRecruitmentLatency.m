
%{ 

This script plots recruitment latencies over the electrodes

%}

dotsize = 70;

% Get the channel locations for the patient
chLocs = Patient.xyChan(:,2:4);
colormin = min(Patient.avgRecruitmentLat);
colormax = max(Patient.avgRecruitmentLat);

% Show all electrode locations, making them unfilled and dark
scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),dotsize,'k')

hold on

for i = 1:size(chLocs,1)
    rLat= Patient.avgRecruitmentLat(i);
    if isnan(rLat) == 0
        scatter3(chLocs(i,1),chLocs(i,2),chLocs(i,3),dotsize,rLat,'filled')
        
    end
    
end

colormap jet
colorbar
