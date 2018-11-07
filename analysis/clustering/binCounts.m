function [t_return,counts] = binCounts(times,window)

total_time = times(end) - times(1);
nbins = ceil(total_time/window);

counts = zeros(nbins,1);
t_return = zeros(nbins,1);

for i = 1:nbins
    bin_times = [times(1) + window*(i-1),min(times(1) + window*i,times(end))];
    t_return(i) = (bin_times(1)+bin_times(2))/2;
    counts(i) = sum(times>=bin_times(1) & times<=bin_times(2));
   
end

end