clear

%% File names
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
%electrodeFile = [electrodeFolder,csvFile];
p1 = genpath(scriptFolder);
addpath(p1);
ptWithSeq = 'ptWithSeq.mat';
finalPt = 'finalPt.mat';

load([resultsFolder,'ptStructs/',finalPt]);

nperm = 1e3;

for i = 1:length(pt)
   for j = 1:length(pt(i).sz)
       
       if isfield(pt(i).sz(j),'data') == 0
          continue 
       end
       if isempty(pt(i).sz(j).data) == 1
           continue
       end
       
       chanData = pt(i).sz(j).data.xyChan;
       sequences = pt(i).sz(j).data.sequences;
       ictal = pt(i).sz(j).data.ictal_idx;
       nictal = sum(ictal);
       
       % Number of sequences
       nseq = pt(i).sz(j).stats.nseqs;
       
       % Get the array of vectors once for real 
       info = getVectors(sequences,nseq,ictal,chanData);
       
       % Get the average vector and get statistics
       avg_ic = mean(info.vectors.ic_vectors);
       avg_interic = mean(info.vectors.interic_vectors);
       angle = acosd(dot(avg_ic,avg_interic)/(norm(avg_ic)*norm(avg_interic)));
       
       % Line length difference
       lldiff = mean(info.lineLength.ic_lineLength) - mean(info.lineLength.interic_lineLength);
        
        %{
        arbitrary_vector = [0,0,1];
        angles_ic = acosd(dot(ic_vectors,repmat(arbitrary_vector,[size(ic_vectors,1),1]))...
            ./(norm(ic_vectors).*norm(repmat(arbitrary_vector,[size(ic_vectors,1),1]))));
        angles_interic = acosd(dot(interic_vectors,repmat(arbitrary_vector,[size(interic_vectors,1),1]))...
            /(norm(interic_vectors)*norm(repmat(arbitrary_vector,[size(interic_vectors,1),1]))));

        std_angles_ic = std(angles_ic);
        std_angles_interic = std(angles_interic);
        
        angle = acosd(dot(avg_ic,avg_interic)/(norm(avg_ic)*norm(avg_interic)));
        %}
        
        all_perm_angles = zeros(nperm,1);
        all_perm_ll = zeros(nperm,1);
        %% Now do a permutation test
        for iperm = 1:nperm
            all_seq = 1:nseq;
            new_ic = randperm(nseq,nictal);
            new_ic_logic = zeros(nseq,1);
            new_ic_logic(new_ic) = 1;
            ictal = new_ic_logic;
            
            info_p = getVectors(sequences,nseq,ictal,chanData);
            avg_ic_p = mean(info_p.vectors.ic_vectors);
            avg_interic_p = mean(info_p.vectors.interic_vectors);
            angle_p = acosd(dot(avg_ic_p,avg_interic_p)/(norm(avg_ic_p)*norm(avg_interic_p)));
            all_perm_angles(iperm) = angle_p;
            
            lldiff_p = mean(info_p.lineLength.ic_lineLength) - mean(info_p.lineLength.interic_lineLength);
            all_perm_ll(iperm) = lldiff_p;
        end
        
        
        
        diffDiff = abs(angle-all_perm_angles);
        [~,I] = min(diffDiff);
        p_vector = 2*min(I/nperm,(1-I/nperm));
        
        diffDiff = abs(lldiff - all_perm_ll);
        [~,I] = min(diffDiff);
        p_ll = 2*min(I/nperm,(1-I/nperm));
        
    end
end