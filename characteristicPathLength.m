function cPL = characteristicPathLength(myNetwork)
% number of nodes
N = numnodes(myNetwork);
% weighted adjacency matrix
% compute shortest path length for every two nodes
shortestPathLength = zeros(N,N);
for i = 1:N
    for j = 1:N
        % shortestpath is smart enough to return zero when the path is
        % from a node to itself! (source = target)
        [~, shortestPathLength(i,j)] = shortestpath(myNetwork,i,j,'method', 'unweighted');
    end
end
cPL = (1/(N*(N-1))) * sum(sum(shortestPathLength));

end