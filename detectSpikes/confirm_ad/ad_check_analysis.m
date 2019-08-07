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

% load patient file
load([structFolder,seq_file])

if isempty(whichPts) == 1
     whichPts = [1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31];
end

for whichPt = whichPts
    ad_rm_spike = nanmean(power(whichPt).alpha.rm_spike./power(whichPt).delta.rm_spike,1);
    ad_keep_spike = nanmean(power(whichPt).alpha.keep_spike./power(whichPt).delta.keep_spike,1);
    
    old_ad = nanmean(old_power(whichPt).alpha./old_power(whichPt).delta,1);
    
    keep_nan = isnan(ad_keep_spike);
    rm_nan = isnan(ad_rm_spike);
    if isequal(keep_nan,rm_nan) == 0
        error('what\n');
    end
    rho = corr(ad_rm_spike(~rm_nan)',ad_keep_spike(~rm_nan)');
    
    figure
    
    plot(ad_rm_spike)
    hold on
    plot(ad_keep_spike)
    plot(old_ad)
    legend('Remove spike','Keep spike','old')
    title(sprintf('%s correlation: %1.4f',pt(whichPt).name,rho))
    pause
    close(gcf)
    
    
end

end