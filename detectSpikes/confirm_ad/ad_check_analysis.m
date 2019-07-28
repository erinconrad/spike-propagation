function ad_check_analysis(whichPts)

%% Load file paths, etc.
[~,~,scriptFolder,resultsFolder,pwfile] = fileLocations;
p1 = genpath(scriptFolder);
addpath(p1);
structFolder = [resultsFolder,'ptStructs/'];

seq_file = 'long_seq.mat';
power_file = 'power_check.mat';


load([structFolder,power_file])

for whichPt = whichPts
    ad_rm_spike = nanmean(power(whichPt).alpha.rm_spike./power(whichPt).delta.rm_spike,1);
    ad_keep_spike = nanmean(power(whichPt).alpha.keep_spike./power(whichPt).delta.keep_spike,1);
    
end

end