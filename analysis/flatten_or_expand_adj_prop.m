function adj_out = flatten_or_expand_adj_prop(adj_in)

%% Detect if it is currently flat or expanded
[m,n] = size(adj_in);
smallest_dim = min([m,n]);
if smallest_dim == 1
    currently_flat = 1;
else
    currently_flat = 0;
end

%% Expand it if it's currently flat
if currently_flat == 1
    
    % Get the number of channels (solving quadratic equation)
    y = length(adj_in);
    nchs = 0.5 + sqrt(0.25+2*y);
    
    % sanity checks
    if nchs-floor(nchs) > 0.01
        error('what\n');
    end
    nchs = round(nchs);
    if nchs*(nchs-1)/2 ~= y
        error('what\n');
    end
    
    % Reconstruct upper triangular matrix
    adj_out = zeros(nchs,nchs);
    count = 0;
    for i = 1:nchs
        for j = 1:i-1
            count = count + 1;
            adj_out(j,i) = adj_in(count);
        end
    end
        
    if count ~= length(adj_in)
        error('what\n');
    end
    
    % Reflect across the diagonal to get full adjacency matrix
    adj_out = adj_out + adj_out';
    
else
    
    %% Flatten it if it's currently expanded
    nchs = size(adj_in,1);
    adj_out = zeros(nchs*(nchs-1)/2,1);
    
    count = 0;
    
    for i = 1:nchs
        for j = 1:i-1
            count = count + 1;
            adj_out(count) = adj_in(j,i);
        end
    end
    
end


end