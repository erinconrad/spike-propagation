function pt = szToIgnore(pt)

for i = 1:length(pt)
    szOnset = zeros(length(pt(i).sz),1);
   for j = 1:length(pt(i).sz)
       pt(i).sz(j).ignoreSz = 0;
       szOnset(j) = pt(i).sz(j).onset;
       
   end
   diffOnset = diff(szOnset);
   for j = 1:length(diffOnset)
      if diffOnset(j) < 12*3600
          pt(i).sz(j+1).ignoreSz = 1;
      end
   end
   
end


end