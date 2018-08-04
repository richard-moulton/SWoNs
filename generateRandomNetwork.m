function randomNetwork = generateRandomNetwork(numNodes,numEdges,randomseed,weightDistribution)
%% generate random adjacency matrix
A = zeros(numNodes,numNodes);
%% create edges
createdEdges = 0;   % keep track of created edges
while createdEdges < numEdges
    addRandomEdge();
end

%% make sure there are no unconnected nodes
% sanity check that no node is unconnected
% (code copied from createNetwork.m)
% added: if we add an edge, remove another one. Do this until all nodes are
% connected but number of edges is preserved (critical for comparison with
% small-world network)

changesMade = 1000; % initially
while changesMade > 0
    changesMade = 0;
    for i=1:numNodes
        % Check for any unconnected rows
       if(sum(A(i,:)) == 0)
          j = i;
          while j == i      % make sure i and j aren't the same
              j = randi([1 numNodes]);
          end
          % found a unique pair of i and j; add edge
          A(i,j) = 1;
          changesMade = changesMade + 1;
          
          %%% need to remove an edge at random to keep # of edges constant
          removeRandomEdge(); 
       end
    end

    for j = 1:numNodes
       % Check for any unconnected columns
       if(sum(A(:,j)) == 0)
          i = j;
          while i == j      % make sure i and j aren't the same
              i = randi([1 numNodes]);
          end
          % found a unique pair; add edge
          A(i,j) = 1;
          changesMade = changesMade + 1;
          
          %%% need to remove an edge at random to keep # of edges constant
          removeRandomEdge();
       end
    end
end

%% give the edges random weights, according to the input distribution
rng(randomseed);
switch weightDistribution
    case 'uniform'
        A_w = A .* rand(numNodes,numNodes);
    case 'gaussian'
        A_w = A .* randn(numNodes,numNodes);
    otherwise 
        disp(['ERROR: UNKNOWN WEIGHT DISTRIBUTION ' weightDistribution]);
end
%% create network object
randomNetwork = digraph(A_w);

function removeRandomEdge()
    % As long as there is an edge to remove in A, removeRandomEdge will
    % randomly select one and remove it.
    
    if sum(sum(A)) == 0
        disp('ERROR: empty adjancency matrix; no edges to remove');
        keyboard;
    end
    
    foundEdge = 0;
    while ~foundEdge
    	newNodes = randperm(numNodes);
        newi = newNodes(1);
        newj = newNodes(2);

        % edge exists; can remove it
        if A(newi,newj)==1
            foundEdge = 1;
        end
      end
      % remove edge that we found
      A(newi,newj) = 0;
end

function addRandomEdge()
    % This function randomly selects two unconnected nodes and inserts a
    % connection into the adjacency matrix.
    i = randi([1 numNodes]);
    j = randi([1 numNodes]);
    if A(i,j)==0 && i ~= j
        A(i,j) = 1;
        createdEdges = createdEdges + 1;
    end
end

end