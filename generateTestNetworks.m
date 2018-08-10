% code is used to generate cell arrays of structs, each containing a
% network object (and associated detail) to pass on to plotting function
% inputs: constant or range of values for each parameter, number of 
% repetitions for each parameter combination, name of .mat file to save this 
% to
% output: cell array of structs, each containing a network object. Also
% saves this to the local directory
% output cell array is N x K x q x numReps; use this to index the array
function networkSet = generateTestNetworks(allN, allK, allq, numReps, saveFileName)
%% allK can be a function of allK; check if this is the case
if isa(allK, 'function_handle')
    % compute N-specific K's
    Ks = allK(allN);
    Kix_range = size(Ks,1);
else    % allK is just a vector of values to be applied to all N's
    Ks = allK;
    Kix_range = length(allK);
end
networkSet = cell(length(allN),length(allK),length(allq),numReps);
% vary N (works also if N is one constant)
for nix = 1:length(allN)
    N = allN(nix);
    disp(['************ N = ' num2str(N)]);
    % vary K (works also if K is one constant)
    for kix = 1:Kix_range
        if isa(allK, 'function_handle')
            % Ks is a matrix with K values dependent on N; extract
            % appropriate K
            K = Ks(kix,nix);
        else
            % Ks is a vector of values to apply to all N
            K = allK(kix);
        end
        disp(['******** K = ' num2str(K)]);
        % vary q (works also if q is one constant)
        for qix = 1:length(allq)
            q = allq(qix);
            disp(['*** q = ' num2str(q)]);
            % repeat numReps times
            for rix = 1:numReps
%                 disp(rix);
                networkSet{nix,kix,qix,rix} = createNetwork(N,K,q,0,0);
            end
        end
    end
end
if isempty(saveFileName) || sum(~isnan(saveFileName))==0
    % create meaningful name to save cell array, if none is provided
    % determine which parameters are constant and which are varying
    
    if length(allN) > 1
        Nstring = 'varyingN_';
    else
        Nstring = ['N' num2str(allN) '_'];
    end
    if length(Ks) > 1
        Kstring = 'varyingK_';
    else
        Kstring = ['K' num2str(K) '_'];
    end
    if length(allq) > 1
        qstring = 'varyingq_';
    else
        qstring = ['q' num2str(q) '_'];
    end
    numRepsString = ['reps' num2str(numReps)];
    
    saveFileName = strcat(Nstring, Kstring, qstring, numRepsString);
end
% save
save(saveFileName, 'networkSet', 'allN', 'allK', 'allq', 'numReps');

end