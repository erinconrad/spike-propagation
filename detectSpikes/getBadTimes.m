function bad = getBadTimes(data)

time_noise = 1000;

bad.empty = [];
bad.noise = zeros(ceil(length(data.values)/time_noise),2);


bad.empty = find(sum(abs(data.values),2)==0);
noise = sqrt(sum((diff(data.values,1,1).^2),2));

for i = 1:length(bad.noise)
  startIdx = (i-1)*time_noise + 1;
  endIdx = min(size(data.values,1)-1,startIdx+time_noise);
  bad.noise(i,1) = (startIdx+endIdx)/2;
  bad.noise(i,2) = mean(noise(startIdx:endIdx));
    
    
end

end