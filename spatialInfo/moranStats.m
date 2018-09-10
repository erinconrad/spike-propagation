% This function takes an array of average recruitment latencies for a group
% of electrodes and a matrix of spatial weights wij and it calculates the
% Moran index, as well as the expected value of the Moran under the null
% distribution, as well as the expected variance (under the asymptotic
% normality assumption)

function MI = moranStats(rl,wij,N)

E = -1/(N-1);


xdmean = rl-nanmean(rl);
W = nansum(wij(:));



op = xdmean' * xdmean;
if sum(size(op)) == 2
    error('Warning, I think you need to take the transpose of rl\n');
end

I = N/W*nansum(nansum(wij.*op))/...
    nansum((xdmean).^2);


%% Alternate way (without matrix multiplication). Confirmed that it gives the same result as above 9/9/18.
%{
topSum = 0;
for i = 1:length(xdmean)
   for j = 1:length(xdmean)
       if isnan(wij(i,j)*xdmean(i)*xdmean(j)) == 0
        topSum = topSum + wij(i,j)*xdmean(i)*xdmean(j);
       end
   end
end
bottomSum = 0;
for i = 1:length(xdmean)
    if isnan((xdmean(i))^2) == 0
        bottomSum = bottomSum + (xdmean(i))^2;
    end
end
I = N/W*topSum/bottomSum;
%}

S1 = 1/2*nansum(nansum((wij+wij').^2));
S2 = nansum((nansum(wij)+nansum(wij')).^2);
S3 = 1/N*nansum((xdmean).^4)/(1/N*nansum((xdmean).^2))^2;
S4 = (N^2-3*N+3)*S1 - N*S2 + 3*W^2;
S5 = (N^2-N)*S1 - 2*N*S2 + 6*W^2;

V = (N*S4-S3*S5)/((N-1)*(N-2)*(N-3)*W^2)-E^2;

Z = (I-E)/sqrt(V);
p = normcdf(-Z);

MI.rl = rl;
MI.I = I;
MI.E = E;
MI.V = V;
MI.Z = Z;
MI.p = p;
MI.wij = wij;
end