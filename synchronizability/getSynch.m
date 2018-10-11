function getSynch(data,Fs)

%% Parameters
freq = [100 150];
window_length = 100;
sample_overlap = 80;
Fs = 1000;

%% Calculate spectral coherence
n_chan = size(data,2);

% Initialize adjacency matrix
A = zeros(n_chan, n_chan);

% Compute all coherences
for i = 1:n_chan
   for j = 1:n_chan
      [out,F] =  mscohere(data(:,i),data(:,j),...
          hamming(window_length),...
          sample_overlap,...
          window_length,...
          Fs);
      
      % Find the indices of the frequency band closest to the desired frequency band
      
       closest_f = intersect(find(F>=freq(1)),find(F<=freq(2)));
       
       % Store coherence in A
       A(i,j) = mean(out(closest_f));
       
   end
    
end

%% Calculate node strength D

% Get the degree vector of the adjacency matrix
dvector = sum(A,1);

% convert it into a diagonal matrix D
D = diag(dvector);

%% Calculate Laplacian L

L = D - A;

%% Compute the eigenspectrum of L(t)
e = eig(L);


%% Calculate synchronizabilty s

% Sort smallest to largest eigenvalue
e = sort(e);

% Compute synchronizability
sync = abs(e(2)/e(end));

end