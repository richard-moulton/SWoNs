function network = createNetwork (N, K, q,displayFlag)
    % Given a number of nodes, N, a degree of connectedness,
    % K, and a rewiring proportion, q, generates and returns
    % a graph object.
    
    %% TO CHANGE BEFORE TESTING:
    % randomize random seed

    % set up random seed for debugging purposes
    randomseed = 5;
    rng(randomseed);
    
    % initialize edges
    A = zeros(N,N);
    
    % connect in ring lattice, with connections between i and it's 2K
    % neighbors (K on either side)
    for i=0:N-1
        for j=1:K
           newi = i+1;
           newj = mod(i+j,N)+1;
           A(newi,newj) = 1;
           newj = mod(i-j,N)+1;
           A(newi,newj) = 1;
        end
    end
    
    if displayFlag
        disp(A);
        % plot
        figure(1);
        subplot(1,3,1);
        plot(digraph(A));
        title('Ring Lattice');
    end
    
    % number of existing edges (assumes all edges are >= 0)
    numEdges = 2*N*K;
    
    % number of edges to rewire
    numRewiredEdges = q*numEdges;
    removedEdges = NaN(numRewiredEdges,2);
    addedEdges = NaN(numRewiredEdges,2);
    
    for e=1:numRewiredEdges        
        % first, find an edge to remove (that hasn't been removed already)
        foundEdge = 0;
        while ~foundEdge
            %% find edge to remove
            % select node at random
            oldi = randi([1 N]);
            % select from its 2K neighbors at random (these will always be
            % edges initially)
            jtemp1 = randi([1 K]);  % how far I'm shifting
            jtemp2 = randi([0 1]); % direction of shift

            if jtemp2==0
                % shift backwards
                oldj = mod(oldi-1 - jtemp1,N)+1;
            else
                % shift forwards
                oldj = mod(oldi-1 + jtemp1,N)+1;
            end

            % make sure this edge hasn't been marked to removed
            if isempty(find(removedEdges(:,1)==oldi & removedEdges(:,2)==oldj,1))
                foundEdge = 1;
            end
        end
        % keep track of removed edges (for debugging)
        removedEdges(e,:) = [oldi oldj];
        % remove edge
        A(oldi, oldj) = 0;
        %%
        
        % select two new nodes to place the edge
        foundEdge = 0;
        while ~foundEdge
            newNodes = randperm(N);
            newi = newNodes(1);
            newj = newNodes(2);
            
            % edge doesn't exist; can put one here
            if A(newi,newj)==0
                foundEdge = 1;
            end
        end
        % keep track of added edges (for debugging)
       addedEdges(e,:) = [newi newj];
        % found a non-existing edge; create it
        A(newi,newj) = 1;
    end
    
    % sanity check that no node is unconnected
    for i=1:N
        % Check for any unconnected rows
       if(sum(A(i,:)) == 0)
          j = i;
          while j == i      % make sure i and j aren't the same
              j = randi([1 N]);
          end
          % found a unique pair of i and j; add edge
          A(i,j) = 1;
       end
    end
    
    for j = 1:N
       % Check for any unconnected columns
       if(sum(A(:,j)) == 0)
          i = j;
          while i == j      % make sure i and j aren't the same
              i = randi([1 N]);
          end
          % found a unique pair; add edge
          A(i,j) = 1;
       end
    end
    
    % for display purposes
    if displayFlag
        disp(removedEdges);
        disp(addedEdges);
        disp(A);
        subplot(1,3,2);
        plot(digraph(A));
        title('Quasi-SWoN');  
    end
    %% assign weights to the edges
    % weights come from a uniform distribution on the open interval (0,1)
    % distribution can be changed, but weights need to remain in the closed
    % interval [0 1] for the math to work, according to Xu 2010
    % keep A as 0's and 1's to indicate presence of edge (we'll use this to
    % create the time delay matrix later); feed in weighted edges directly
    % into the digraph function to create our network
    
    % create network object
    newNetwork = digraph(rand(N,N) .* A);
    
    %% measuring small-world-ness
    % small-world-ness S = 
    % to measure small-world-ness, we will use measures detailed in Xu 2010
    % these require that we calculate the shortest path length between each
    % two nodes, and the clustering coefficient for our network, and a
    % random network
    [charPathLength, clusterCoeff] = networkStats(newNetwork);
    % do the same for a randomly generated network
    randomNetwork = generateRandomNetwork(N,numedges(newNetwork),randomseed,'uniform');
    subplot(1,3,3);
    plot(randomNetwork);
    title('Random');
    [charPathLength_random, clusterCoeff_random] = networkStats(randomNetwork);
    
    % ratio of characteristic path lengths
    lambda = charPathLength / charPathLength_random;
    % ratio of clustering coefficients
    gamma = clusterCoeff / clusterCoeff_random;
    % small-world-ness!
    smallWorldMeasure = gamma / lambda; 
    
    if displayFlag
        disp([smallWorldMeasure gamma lambda]);    
    end
end