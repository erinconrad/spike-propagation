function P = countGridsStripsDepths(P)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile] = fileLocations;
filename = 'pWithElectrodes.mat';

all_strips = [];
all_grids = [];
all_depths = [];

for i = 1:length(P)
    
   n_strips = 0;
   strip_idx = [];
   n_grids = 0;
   grid_idx = [];
   n_depths = 0;
   depth_idx = [];
   n_other = 0;
   other_idx = [];
   
   if isempty(P(i).electrodeData) == 0 && isfield(P(i).electrodeData,'electrodes') == 1
     
   
   for j = 1:length(P(i).electrodeData.electrodes)
       if strcmp(P(i).electrodeData.electrodes(j).type,'S') == 1
           n_strips = n_strips + 1;
           strip_idx = [strip_idx,j];
       elseif strcmp(P(i).electrodeData.electrodes(j).type,'G') == 1
           n_grids = n_grids + 1;
           grid_idx = [grid_idx,j];
       elseif strcmp(P(i).electrodeData.electrodes(j).type,'D') == 1
           n_depths = n_depths + 1;
           depth_idx = [depth_idx,j];
           
           fprintf('%s has depth electrode %s\n',P(i).name,P(i).electrodeData.electrodes(j).name);
       else
           n_other = n_other + 1;
           other_idx = [other_idx,j];
           fprintf('Surprising electrode type for patient %d and electrode %d\n',i,j);
       end
       
   end
   
   end
   
   P(i).electrodeData.n_strips = n_strips;
   P(i).electrodeData.strip_idx = strip_idx;
   
   P(i).electrodeData.n_grids = n_grids;
   P(i).electrodeData.grid_idx = grid_idx;
   
   P(i).electrodeData.n_depths = n_depths;
   P(i).electrodeData.depth_idx = depth_idx;
   
   all_strips = [all_strips;n_strips];
   all_grids = [all_grids;n_grids];
   all_depths = [all_depths;n_depths];
    
end

all_strips
all_grids
all_depths

save([resultsFolder,filename])

end