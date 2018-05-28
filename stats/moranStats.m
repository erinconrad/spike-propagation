% This function takes an array of average recruitment latencies for a group
% of electrodes and a matrix of spatial weights wij and it calculates the
% Moran index, as well as the expected value of the Moran under the null
% distribution, as well as the expected variance (under the asymptotic
% normality assumption)

function MI = moranStats(rl,wij,N)

E = -1/(N-1);


xdmean = rl-nanmean(rl);
W = sum(wij(:));

op = xdmean' * xdmean;
I = N/W*nansum(nansum(wij.*op))/...
    nansum((xdmean).^2);

S1 = 1/2*sum(sum((wij+wij').^2));
S2 = sum((sum(wij)+sum(wij')).^2);
S3 = 1/N*sum((xdmean).^4)/(1/N*sum((xdmean).^2))^2;
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