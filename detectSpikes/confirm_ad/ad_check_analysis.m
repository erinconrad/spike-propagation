function ad_check_analysis(whichPts)

%{
This takes the structure power_check.mat, generated from ad_remove_spikes,
to compare the alpha-delta ratio when I keep spikes to that when I remove
spikes. It also compares both to my original calculation of the alpha/delta
ratio. 
%}

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
    % The patients I studied in the paper
     whichPts = [1,4,6,7,8,9,12,14,15,16,17,18,19,20,22,24,25,27,30,31];
end

for whichPt = whichPts
    ad_rm_spike = nanmean(power(whichPt).alpha.rm_spike./power(whichPt).delta.rm_spike,1);
    ad_keep_spike = nanmean(power(whichPt).alpha.keep_spike./power(whichPt).delta.keep_spike,1);
    
    % Get my old alpha-delta ratio (which I calculated in a slightly
    % different way) to make sure that it isn't different from my new
    % alpha-delta ratios (both the ones where I keep and where I remove
    % spikes)
    old_ad = nanmean(old_power(whichPt).alpha./old_power(whichPt).delta,1);
    
    keep_nan = isnan(ad_keep_spike);
    rm_nan = isnan(ad_rm_spike);
    if isequal(keep_nan,rm_nan) == 0
        error('what\n');
    end
    rho = corr(ad_rm_spike(~rm_nan)',ad_keep_spike(~rm_nan)');
    rho_old = corr(ad_rm_spike(~rm_nan)',old_ad(~rm_nan)');
    
    figure
    
    plot(ad_rm_spike)
    hold on
    plot(ad_keep_spike)
    plot(old_ad)
    legend('Remove spike','Keep spike','old')
    title(sprintf('%s rm/keep correlation: %1.4f, rm/old correlation: %1.4f',pt(whichPt).name,rho,rho_old))
    pause
    close(gcf)
    
    
end

end