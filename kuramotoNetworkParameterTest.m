%% Test all combinations of network and Kuramoto oscillator parameters 
% N, K, q, and lam may be set to single values or ranges, and the script
% will test all resulting combinations. Note that this will result in an
% infeasible computation time if too many ranges are given.
% TODO: variable omega and theta0
clear;
multiWaitbar('CloseAll');

% flags etc.
displayFlag = false;  % false = suppress plots
saveNetworkFlag = false;
overallProgressFmt = 'Overall progress (current N=%d)';
netProgressLabel = 'Current network progress';

% network parameters 
N = [20];  %[2 5 10 25 50 100 250];  % number of nodes
Kprop = [0.4];  % edge density parameter; 0 -> K=0 and 1 -> K=N/2
q = [0.1];  % rewiring parameter
netParamCombs = num2cell(allcomb(N, Kprop, q));  % get combinations

% oscillator parameters (maybe more of them in future)
lam = linspace(0, 2, 50); %[0, logspace(-7, 2, 5)];
oscParamCombs = num2cell(allcomb(lam));  % get combinations

% simulation parameters
h = 1e-1;  % numerical integration step
steps = 500;

% pre-allocate arrays
% theta changes in size depending on network...
theta = cell(size(netParamCombs, 1), size(oscParamCombs, 1));
% z is a scalar per step, no matter the network parameters
% r and psi can be recovered later, no need to store everything
z = zeros(size(netParamCombs, 1), size(oscParamCombs, 1), steps);

% hide any internally generated plots and setup progress bars
set(0, 'defaultfigurevisible', 'off');
overallProgressLabel = sprintf(overallProgressFmt, '');
multiWaitbar(overallProgressLabel, 0);
multiWaitbar(netProgressLabel, 0);

% iterate over all combinations of network parameters
for j = 1:size(netParamCombs, 1)
    [N, Kprop, q] = deal(netParamCombs{j, :});

    % initial/constant conditions
    omega = rand(1, N);  %initialize nodes with random intrinsic frequency
    theta0 = 2 * pi * rand(1, N);

    % create network
    K = floor(Kprop * N / 2);
    % TODO: don't repeat if K or q are (effectively) 0 more than once
    %if K == 0 && find(netParamCombs(1:j-1))
    %    continue
    %end
    net = createNetwork(N, K, q, displayFlag, saveNetworkFlag);

    % update overall waitbar
    tmpLabel = sprintf(overallProgressFmt, N);
    if ~strcmp(tmpLabel, overallProgressLabel)
        f = multiWaitbar(overallProgressLabel, (j - 1)/size(netParamCombs, 1), 'Relabel', tmpLabel);
        overallProgressLabel = tmpLabel;
    end

    % iterate over all oscillator parameter combinations
    for i = 1:size(oscParamCombs, 1)
        [lam] = deal(oscParamCombs{i,:});
        multiWaitbar(netProgressLabel, i/length(oscParamCombs));
        [theta{j, i}, z(j, i, :)] = kuramNetwork(net, lam, omega, theta0, steps, h);
    end
end

% close waitbar and re-enable figure visibility
multiWaitbar('CloseAll');
set(0,'defaultfigurevisible','on');




