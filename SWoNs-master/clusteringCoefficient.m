function C = clusteringCoefficient(myNetwork)
% number of nodes
N = numnodes(myNetwork);
% weighted adjacency table
weightedEdges = table2array(myNetwork.Edges);
% convert to matrix
A_w = zeros(N,N);
for e = 1:size(weightedEdges,1)
    A_w(weightedEdges(e,1), weightedEdges(e,2)) = weightedEdges(e,3);
end
        
% compute characteristic path length per node
c = NaN(N,1);
for i = 1:N
    num = 0;
    denom1 = 0;
    denom2 = 0;
    denom3 = 0;
    for j = 1:N
        for k = 1:N
            num = num + A_w(i,j) * A_w(j,k) * A_w(k,i);
        end
        denom1 = denom1 + A_w(i,j);
        denom2 = denom2 + A_w(j,i);
        denom3 = denom3 + A_w(i,j) * A_w(j,i);
    end
    c(i) = num / (denom1*denom2 - denom3);
end
% average network clustering coefficient
C = sum(c) / N;
end