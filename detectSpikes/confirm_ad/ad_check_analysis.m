function ad_check_analysis(whichPts)

%% Load file paths, etc.
[~,~,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];

seq_file = 'long_seq.mat';
power_file = 'power_check.mat';
old_power_file = 'power2.mat';

power=load([structFolder,power_file]);
power=power.power;

old_power = load([structFolder,old_power_file]);
old_power = old_power.power;

for whichPt = whichPts
    ad_rm_spike = nanmean(power(whichPt).alpha.rm_spike./power(whichPt).delta.rm_spike,1);
    ad_keep_spike = nanmean(power(whichPt).alpha.keep_spike./power(whichPt).delta.keep_spike,1);
    
    old_ad = nanmean(old_power(whichPt).alpha./old_power(whichPt).delta,1);
    
    rho = corr(ad_rm_spike',ad_keep_spike');
    
    figure
    
    plot(ad_rm_spike)
    hold on
    plot(ad_keep_spike)
    plot(old_ad)
    legend('Remove spike','Keep spike','old')
    title(sprintf('Correlation: %1.2f',rho))
    pause
    close(gcf)
    
    
end

end