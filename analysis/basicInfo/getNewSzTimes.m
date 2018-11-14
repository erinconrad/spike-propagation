function pt = getNewSzTimes(pt,ptNew)

for i = 1:length(pt)
   for inew = 1:length(ptNew)
      if strcmp(pt(i).name,ptNew(inew).name) == 1
         pt(i).newSzTimes = [];
         for j = 1:length(ptNew(inew).sz)
            pt(i).newSzTimes = [pt(i).newSzTimes;...
                ptNew(inew).sz(j).onset ptNew(inew).sz(j).offset];
               
          
         end
      end
       
       
       
   end
    
    
end