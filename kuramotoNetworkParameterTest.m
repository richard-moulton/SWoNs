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
netGenReps = 10; % number of times to generate network with same parameters (e.g. for non-zero q/non-deterministic network structure)
netReps = 1;  % number of times to run same network (e.g. to average over noise)
N_vec = [20];  %[2 5 10 25 50 100 250];  % number of nodes
q = [0]; % 0.05 0.1 0.15]; %linspace(0,1,50)  % rewiring parameter
K_fcn = @(N) linspace(0, 1, floor(N/2) + 1) * floor(N/2);
netParamCombs = [];
for i = 1:length(N_vec)
   K_vec = K_fcn(N_vec(i));
   netParamCombs = [netParamCombs, allcomb(N_vec(i), K_vec, q)]; 
end
netParamCombs = num2cell(netParamCombs);

% oscillator parameters (maybe more of them in future)
%fac = gamma(2:0.1:14 + 1);
lam = linspace(0, 1, 100); %[0 1./fac(end:-1:1)]; %logspace(-8, 1, 100)]; %  
oscParamCombs = num2cell(allcomb(lam));  % get combinations

% simulation parameters
h = 1e-2;  % numerical integration step
steps = 300;

% pre-allocate arrays
% theta changes in size depending on network...
nets = cell(size(netParamCombs, 1));
theta = cell(size(netParamCombs, 1), netGenReps, netReps, size(oscParamCombs, 1));
% z is a scalar per step, no matter the network parameters
% r and psi can be recovered later, no need to store everything
z = zeros(size(netParamCombs, 1), netGenReps, netReps, size(oscParamCombs, 1), steps);

% hide any internally generated plots and setup progress bars
set(0, 'defaultfigurevisible', 'off');
overallProgressLabel = sprintf(overallProgressFmt, '');
multiWaitbar(overallProgressLabel, 0);
multiWaitbar(netProgressLabel, 0);

% iterate over all combinations of network parameters
for j = 1:size(netParamCombs, 1)
    [N, K, q] = deal(netParamCombs{j, :});

    % initial/constant conditions
    omega = rand(1, N);  %initialize nodes with random intrinsic frequency
    theta0 = 2 * pi * rand(1, N);

    % update overall waitbar
    f = multiWaitbar(overallProgressLabel, (j - 1)/size(netParamCombs, 1));
    tmpLabel = sprintf(overallProgressFmt, N);
    if ~strcmp(tmpLabel, overallProgressLabel)
        f = multiWaitbar(overallProgressLabel, 'Relabel', tmpLabel);
        overallProgressLabel = tmpLabel;
    end

    for k = 1:netGenReps  % number of networks to (stochastically) generate
    % iterate over all oscillator parameter combinations
        nets{j, k} = createNetwork(N, K, q, displayFlag, saveNetworkFlag);
        for l = 1:netReps  % number of repetitions per (noisy) network
            % TODO: move netReps inside kuramNetwork/dimension in theta
            % array?
            for i = 1:size(oscParamCombs, 1)
                [lam] = deal(oscParamCombs{i,:});
                multiWaitbar(netProgressLabel, i/length(oscParamCombs));
                [theta{j, k, l, i}, z(j, k, l, i, :)] = kuramNetwork(nets{j, k}, lam, omega, theta0, steps, h);
            end
        end
    end
end

% close waitbar and re-enable figure visibility
multiWaitbar('CloseAll');
set(0,'defaultfigurevisible','on');




