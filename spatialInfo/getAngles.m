function angles = getAngles(sequences,electrodeData)

nseq = size(sequences,2);
angles = zeros(nseq,1);

%% Loop through sequences
for i = 1:nseq
   curr_seq = sequences(:,i); 
   [B,I] = sort(curr_seq);
   C = I(isnan(B) == 0);
   
   nspikes = length(C);
   
   % get locations of the spikes
   locs = electrodeData.locs(C,2:4);
   
   
   % if odd number of spikes, ignore the middle spike
   if mod(nspikes,2) == 1
       early_idx = 1:(nspikes-1)/2;
       late_idx = (nspikes-1)/2+2:nspikes;
   else
       early_idx = 1:nspikes/2;
       late_idx = nspikes/2+1:nspikes;
   end
    
   % get mean coordinates of early channels and late channels
   early_mean = mean(locs(early_idx,:),1);
   late_mean = mean(locs(late_idx,:),1);
   
   % get vector between them 
   vec = late_mean - early_mean;
   
   ref_vector_both = electrodeData.ref_vector;
   ref_vector = ref_vector_both(2,:) - ref_vector_both(1,:);
   
   % get angle between this vector and reference vector
   angle = acosd(dot(vec,ref_vector)/norm(vec)/norm(ref_vector));
   
   angles(i) = angle;
   
   
   % sample plot
   %{
   chLocs = electrodeData.locs(:,2:4);
   figure
   scatter3(chLocs(:,1),chLocs(:,2),chLocs(:,3),...
        30,'b','filled')
    hold on
   scatter3(early_mean(1),early_mean(2),early_mean(3),30,'g','filled')
   scatter3(late_mean(1),late_mean(2),late_mean(3),30,'r','filled')
   
   scatter3(ref_vector_both(1,1),ref_vector_both(1,2),...
       ref_vector_both(1,3),30,'g','filled')
   scatter3(ref_vector_both(2,1),ref_vector_both(2,2),...
       ref_vector_both(2,3),30,'r','filled')
   
   %}
   
end

end