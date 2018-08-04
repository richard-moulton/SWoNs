function network = createNetwork (N, K, type,displayFlag)
    % Given a number of nodes, N, a degree of connectedness,
    % K, and a desired type (random, small-world, scale
    % free), generates and returns a graph object.
    
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
        subplot(1,2,1);
        plot(digraph(A));
        title('Ring Lattice');
    end
    
    % number of existing edges (assumes all edges are >= 0)
    numEdges = 2*N*K;
    
    % rewiring probability (between 0 and 1)
    q = 0.2;
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
    
    A(1,:) = 0;
    A(:,5) = 0;
    
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
        subplot(1,2,2);
        plot(digraph(A));
        title('Quasi-SWoN');  
    end
    
    % create network object
    network = digraph(A);
    keyboard;
end