% function generates various plots characterizing the relationship between
% different network parameters
% inputs: cell array of networks to pass on to the oscillator, label of
% plot to generate
% outputs: saves plots to the current directory with hopefully meaningful
% names
% set lambda values to test on in code
% note: user's responsibility to enter the appropriate network set! code
% will double-check for this
%% to run on Matlab 2018; either change 'kuramNetwork2015' to 'kuramNetwork', or download 'kuramNetwork2015' (it's forward compatible)
function generateParameterPlots2015(networkSetFile, plotID)
% allLambdas = linspace(0,1,5);
allLambdas = linspace(0,0.5,100);
h = 1e-2;  % for kuramoto oscillator: numerical integration step
syncThreshold = 0.99;   % threshold beyond which to define synchrony
timeSteps = 300;    % for kuramoto oscillator
% omega_std = 1;      % standard deviation for the distribution of omega (intrinstic frequencies)
% theta0_std = 1;     % standard deviation for the distribution of initial theta (phases)

%% load provided network set file
networkSetElements = load(networkSetFile);
%% extract elements and name them
networkSet = networkSetElements.networkSet;
allN = networkSetElements.allN;
allK = networkSetElements.allK;
allq = networkSetElements.allq;
numReps = networkSetElements.numReps;

%% determine which plot we're making
switch plotID
    case 1  % heatmap of R values for 
        plot1;
    case 2  % plot time to synchrony vs. time for different values of K
        plot2;
    otherwise
        disp(['UNKNOWN PLOT ID ' num2str(plotID)]);
        keyboard;
end

    function plot1
        %% make sure network set is appropriate for this plot
        % should now have the range of N, K and q values and the number of
        % repetitions for each
        % N, K and q should be fixed for these plots; numReps can vary
        
        if length(allN) ~= 1 || length(allK) ~= 1 || length(allq) ~= 1
            disp('ERROR: inappropriate network set');
            keyboard;
        else
            N = allN;   % only need N
%             K = allK;
%             q = allq;
        end
        
        % initialize intrinsic frequencies
        omega = rand(1, N);     % uniform distribution
        % initialize phases
        theta0 = 2 * pi * rand(1, N);   % uniform distribution
        
        % compute synchrony values at each time step from 1 to timeSteps,
        % for each network repetition
        allRValues = NaN(length(allLambdas),timeSteps,numReps);
        % mean and standard error of the mean over repetitions (if there are any)
        mean_Rvalues = NaN(length(allLambdas), timeSteps);
        sem_Rvalues = NaN(length(allLambdas), timeSteps);
        for lix = 1:length(allLambdas)
            lambda = allLambdas(lix);
            for r = 1:numReps
                % initialize intrinsic frequencies
                omega = rand(1, N);     % uniform distribution
                % initialize phases
                theta0 = 2 * pi * rand(1, N);   % uniform distribution
                [~, zValues] = kuramNetwork2015(networkSet{1,1,1,r}.networkObject, lambda, omega, theta0, timeSteps, h);
                allRValues(lix,:,r) = computeR(zValues);
            end
            mean_Rvalues(lix,:) = mean(allRValues(lix,:,:),3);
            sem_Rvalues(lix,:) = std(allRValues(lix,:,:),[],3) / sqrt(numReps);
        end
        % plot!
        imagesc(mean_Rvalues);
        colorbar;
    end
end

function r = computeR(z)
    r = abs(z);
end
