function [glme,lme] =  statsSeqFreq(pt,whichPts,window)

%% Remaining problems
% I have no way to remove seizures from the analysis if there are large
% chunks of time in the 6 hours prior to the seizure where data is missing

% Should think about whether the MI would be normally distributed within a
% window

% Seems to be very sensitive to the choice of window :(


%% Remove EKG artifact and depth electrodes
rmEKG = 1;
prox = 0.02; %20 ms
rmDepth = 1;
rmType = 'D';

%% Parameters
nwindows = 6*3600/window;

% The number of the seizure index
szIndex = 0;

% Number of sequences in each window
num_sequences = [];

% An mx1 array of all patient numbers for each window where m is the number
% of windows
pat = [];

% mx1 array of all seizure index numbers
sz = [];

% mx1 array of all window numbers
chunk = [];

% mx1 array of MI
MI = [];

%% Loop through the patients
for i = whichPts
    
   %% Get wij
    xyChan = pt(i).electrodeData.locs;
    dmin = pt(i).dmin;
    wij = getwij(xyChan,dmin);
    nchs = length(pt(i).channels);
    
   for j = 1:length(pt(i).sz)
       
       % Get sequence
       seq = pt(i).sz(j).seq_matrix;
       
       % Get first spikes in each sequence
       first_spikes = min(seq,[],1);
       
       % Skip if too close to the beginning of the file
        if pt(i).sz(j).onset < window*nwindows
            continue
        end
        
        % Skip if too close to the previous seizure
        if j > 2 && pt(i).sz(j).onset - pt(i).sz(j-1).onset < window*nwindows
            continue
        end
       
        % Skip this seizure if the first spike sequence is within 5 hours of
        % the seizure (happens when there are periods of disconnection, e.g.
        % HUP070)
        if pt(i).sz(j).onset - first_spikes(1) < window*(nwindows-1)
            continue
        end
       
        % if haven't advanced in the loop by now, we're keeping the seizure
        szIndex = szIndex + 1;
        
        %% Now get window data
        startTime = pt(i).sz(j).onset - window*nwindows;
        
        % Loop through all the windows
        for tt = 1:nwindows
            
            % Get the times for the windows
            times = [(tt-1)*window + startTime,tt*window + startTime];
           
            
            % Get the number of sequences in the window
            num_sequences = [num_sequences;...
                sum(first_spikes >= times(1) & first_spikes <= times(2))];
            
            % Get the patient id for that window
            pat = [pat;i];
            
            % get the seizure id for that window
            sz = [sz;szIndex];
            
            % get the chunk (or window) id for that window
            chunk = [chunk;tt];
            
            %% Also get MI for the window
            % Get the appropriate sequences for the window
            correct_seqs = seq(:,first_spikes >= times(1) & first_spikes <= times(2));
            
            % Get the latency for each sequence
            latency_all_seq = correct_seqs - min(correct_seqs,[],1);
            
            % Take the average latency for the channel over all sequences
            mean_latency = nanmean(latency_all_seq,2);
            rl = mean_latency';
            
            % Get the moran index
            MIstruct= moranStats(rl,wij,nchs);
            MI = [MI;MIstruct.I];
            
        end
        
        
       
   end
    
end

%% Poisson regression to combine data for all spikes and all seizures
% http://math.bu.edu/people/mak/samsi/SAMSI_GLM_Example.pdf

% This method assumes that spikes from all
% patients and all seizures come from a population with the same mean. This
% is incorrect and will likely cause problems as I include more patients
% and more seizures. However, for small numbers of seizures it appears to
% produce the same result as a generalized linear mixed model.
[b,dev,stats] = glmfit(chunk,num_sequences,'poisson');

%% Plot the data against the model

lambda = exp(b(1) + b(2)*chunk);
figure
scatter(chunk,num_sequences,'b')
hold on
scatter(chunk,lambda,'r')
xlabel('which window')
ylabel('sequence count')
legend('Data','GLM');


%% Make a table for the purpose of doing a generalized linear mixed model for seq freq
T = table(num_sequences,chunk,sz,pat,'VariableNames',{'num_sequences','chunk','sz','pat'});
glme = fitglme(T,'num_sequences ~ 1 + chunk + (1|sz) + (1|pat)',...
    'Distribution','Poisson','Link','log','FitMethod','Laplace',...
    'DummyVarCoding','effects');

%% Make a table for linear mixed model for MI
T2 = table(MI,chunk,sz,pat,'VariableNames',{'MI','chunk','sz','pat'});
lme = fitlme(T2,'MI ~ 1 + chunk + (1|sz) + (1|pat)');

% do a basic linear regression for plotting purposes
[b_MI,~,~,~,stats_MI] = regress(MI,[ones(length(chunk),1),chunk]);
figure
scatter(chunk,MI,'b')
hold on
scatter(chunk,b_MI(1) + b_MI(2)*chunk,'r');
xlabel('which window')
ylabel('MI')
legend('Data','linear regression');

end