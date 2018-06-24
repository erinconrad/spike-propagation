function [ seq ] = spt_initseq( xlocs,ylocs,zlocs,tstamps,idx )
% This function initializes the reference sequence

    X   = xlocs(:,idx);   X(X==0) = [];
    Y   = ylocs(:,idx);   Y(Y==0) = [];
    Z   = zlocs(:,idx);   Z(Z==0) = [];
    T   = tstamps(:,idx); T = T(1:length(X));
    seq = [X,Y,Z,T];

end

    
