%% removeChs

% This function removes spikes that occur in a channel of a specific type,
% so, for instance, we could remove all spikes that occur on depth
% electrodes

function new_gdf = removeChs(gdf,electrodeData,type)

keep_idx= true(size(gdf,1),1);

for i = 1:size(gdf,1)
   ch = round(gdf(i,1));
   ch_type = electrodeData.electrodes(ch).type;
   if ismember(ch_type,type) == 1
       keep_idx(i) = false;
   end

end

new_gdf = gdf(keep_idx,:);


%{
for i = 1:size(new_gdf)
   ch = round(new_gdf(i,1));
   ch_type = electrodeData.electrodes(ch).type;
   if ismember(ch_type,type) == 1
       fprintf('Spike %d electrode type is %s\n',i,ch_type);
   end
    
end
%}

end