function compareStructs(pt1,pt2)

for i = 1:length(pt1)
    if strcmp(pt1(i).name,pt2(i).name) == 0
        fprintf('Warning, look at names for %s\n',pt1(i).name);
    end
    
    if length(pt1(i).sz) ~= length(pt2(i).sz)
        fprintf('Warning, look at number of seizures for %s\n',pt1(i).name);
        continue
    end
    
   for j = 1:length(pt1(i).sz)
      if pt1(i).sz(j).onset ~= pt2(i).sz(j).onset 
        fprintf('Warning, look at onset times for %s seizure %d\n',pt1(i).name,j);
      end
       
      if isequal(pt1(i).sz(j).electrodes, pt2(i).sz(j).electrodes) == 0
          fprintf('Warning, look at onset chs for %s seizure %d\n',pt1(i).name,j);
      end
    
   end
end