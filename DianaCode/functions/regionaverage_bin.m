function [xunique,valunique] = regionaverage_bin(xin,valin,n)
% this function projects data with repeated x positions onto a unique x grid, 
% averaging values at the same x positions

[Y,E] = discretize(xin,n); % Y gives the bin index for each data point
freqs = accumarray(Y',1); % how many entries in each bin 
% [xunique,ia,ic] = unique(xin); %output is in sorted order
% freqs = accumarray(ic,1); % how often does each value occur
% nunique = length(xunique);
 nvals = length(valin);

% matrix which associates with positions in xunique the values of valin (with frequencies as weights)
%Bconnect = sparse(ic,1:nvals,ones(size(ic)),nunique,nvals,nvals);
Bconnect = sparse(Y,1:nvals,ones(size(Y)),n,nvals,nvals);
valunique = ((Bconnect*valin')./freqs)'; % this takes the average value

xunique = ((Bconnect*xin')./freqs)';%0.5*(E(1:end-1)+E(2:end));

% figure(300)
% scatter(xin,valin,'blue','o','MarkerEdgeAlpha',0.2);
% hold on
% plot(xunique,valunique,'b.-');

end