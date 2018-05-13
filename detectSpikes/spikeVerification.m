%{
This function loops through desired patients and loops through various
tmuls and absthresh's and runs portVisualizeSpikes to see how the spike
detector works for that time period

%}


function spikeVerification(P)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;

startTime = 632555.81;
thresholds_to_try = [5 13 18];
absthresh_to_try = [50 300 600];
times_to_try = [3];
soz = [0];
chIds = [10 20 29 30 40];

for i = 1%:length(P)
    outputFolder = [resultsFolder,'spike verification/',P(i).name,'/'];

    if isempty(dir([outputFolder,'*.png'])) == 0
        fprintf('Patient %s presumably done, skipping\n',P(i).name);
        %continue
    end
    
    if isempty(P(i).ieeg_name) == 1 || isempty(P(i).electrode_labels) == 1
        fprintf('Missing ieeg_name or electrode labels for patient %s, skipping\n',P(i).name);
        continue
    end
    
    for s = soz
    
        nchs = size(P(i).electrodeData.locs,1);
        %chIds = randsample(nchs,10);
    for j = times_to_try
       
        for k = thresholds_to_try
            
            for m = absthresh_to_try
            
                portVisualizeSpikes(P,i,1,j,s,chIds,k,m,startTime);
                
            end
            
            
            
        end
        
    end
    
    end
    
    
end



end