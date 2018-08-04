% returns characteristic path length and clustering coefficient for a given
% network object
function [charPathLength, clusterCoeff] = networkStats(myNetwork)
charPathLength = characteristicPathLength(myNetwork);
clusterCoeff = clusteringCoefficient(myNetwork);
end