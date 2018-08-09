function [thetas, zs, netParamCombs, oscParamCombs] = testKuramotoNetParams(N_vec, q_vec, K_fcn, lam_vec, netReps, netGenReps, steps, h)
%testKuramotoNetParams Test combinations of network and Kuramoto oscillator parameters. 
% N_vec, q_vec, and lam_vec may be given as single values or ranges. K_fcn 
% is a function that generates one or more K values to test, and may take
% N as its sole argument.
%
% Note that this will result in an unfeasible computation time if too many 
% large parameter ranges are given.
%
% TODO: variable omega and theta0

clear;
multiWaitbar('CloseAll');

% flags etc.
displayFlag = false;  % false = suppress plots
saveNetworkFlag = false;
overallProgressFmt = 'Overall progress (current N=%d)';
netProgressLabel = 'Current network progress';

% network parameters 
netParamCombs = [];
for i = 1:length(N_vec)
   K_vec = K_fcn(N_vec(i));
   netParamCombs = [netParamCombs, allcomb(N_vec(i), K_vec, q_vec)]; 
end
netParamCombs = num2cell(netParamCombs);

% oscillator parameters (maybe more of them in future)
oscParamCombs = num2cell(allcomb(lam_vec));  % get combinations

% pre-allocate arrays
% theta changes in size depending on network...
nets = cell(size(netParamCombs, 1));
thetas = cell(size(netParamCombs, 1), netGenReps, netReps, size(oscParamCombs, 1));
% z is a scalar per step, no matter the network parameters
% r and psi can be recovered later, no need to store everything
zs = zeros(size(netParamCombs, 1), netGenReps, netReps, size(oscParamCombs, 1), steps);

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
                [thetas{j, k, l, i}, zs(j, k, l, i, :)] = kuramNetwork(nets{j, k}, lam, omega, theta0, steps, h);
            end
        end
    end
end

% close waitbar and re-enable figure visibility
multiWaitbar('CloseAll');
set(0,'defaultfigurevisible','on');

end



