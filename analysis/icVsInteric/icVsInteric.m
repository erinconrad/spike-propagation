function icVsInteric(ptIc,ptInter,whichPts)

for whichPt = whichPts

%% Get sequences

seq_inter = ptInter(whichPt).seq_matrix;
seq_ic = ptIc(whichPt).seq_matrix;


%% Remove interictal sequences that are during the seizure time

%% Remove sequences with too many ties???

%% Group all sequences
all_seq = [seq_inter,seq_ic];

%% Designate ic as a way to keep track of what is ictal as I will delete things later
ic = [zeros(1,size(seq_inter,2)),ones(1,size(seq_inter,1))];

%% Get lead channels

%% Use elbow method to get optimal cluster number

%% Cluster according to location of lead channel

%% Remove clusters that appear to be mostly artifact (because otherwise might just say that artifact more likely to occur during or not during the seizure)

%% Do chi2 to determine whether the cluster distribution is different between the ictal and interictal spikes
[tbl,chi2,p,labels] = crosstab([idx_inter,idx_ic],ic);

%% Do bootstrap to validate chi2
% Randomize ictal versus interictal designation, and see whether ?????
nb = 1e3;
for ib = 1:nb

end


end


end
