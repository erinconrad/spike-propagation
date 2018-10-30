function getDataDuration(pt)

for i = 1:length(pt)
   fprintf('%s\n',pt(i).name); 
   if isempty(pt(i).sz) == 1
       continue
   end
   firstSzOnset = pt(i).sz(1).onset;
   lastSzOffset = pt(i).sz(end).offset;
   fprintf('Total time between seizures: %1.1f hours\n',...
       (lastSzOffset-firstSzOnset)/3600);
end


end