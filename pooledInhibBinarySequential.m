function [thetas, r, z, lams, decision] = pooledInhibBinarySequential(nets, trueSeq, thres, alpha, beta, delta, phi, forgetPenalty, eps, omegas, thetas0, lams0, maxsteps, h)
% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.

% number of trials in sequence
nSeq = length(trueSeq);

% number of networks
nNets = size(nets, 2);

% prior 
lamsPrior = lams0;

% network parameters
for b = 1:nNets
    N(b) = nets{b}.numnodes;
    E{b} = weightedA(nets{b});
    
    % pre-allocate arrays and set initial conditions
    thetas{b} = zeros(N(b), nSeq, maxsteps);
end

% initialize decision matrix
decision = zeros(nNets, nSeq, 2);

exc = [2 1];

for s = 1:nSeq
    lams(:, s, 1) = lamsPrior;
    for b = 1:nNets
        thetas{b}(:, s, 1) = thetas0{b};
        z(b, s, 1) = sum(exp(1i*thetas{b}(:, s, 1))) / N(b); 
        r(b, s, 1) = abs(z(b, s, 1)); 
    end
    decisionTrial(s);
    thisDecision = squeeze(decision(:, s, :));
    if sum(logical(thisDecision(:, 1))) == 1
        decisionStep = max(thisDecision(:, 2));
        whichDecision = find(thisDecision(:, 1));
        correct = 2*(whichDecision == trueSeq(s)) - 1;
        for b = 1:nNets
            decay = phi;
            if b == whichDecision && correct < 0
                decay = forgetPenalty * decay;
            end
            lamsPrior(b) = decay*lamsPrior(b) + correct*delta*log(r(b, s, decisionStep) / r(exc(b), s, decisionStep));
        end
    end
end

function [] = decisionTrial(s)
    % separated into function to preserve any(decision) return 
    for iter = 1:maxsteps-1
        for b = 1:nNets
            % pairwise inter-node phase differences
            thetaPairwiseDiffs = thetas{b}(:, s, iter) - thetas{b}(:, s, iter)'; 

            % numerical integration step
            dw = randn * eps;
            lams(b, s, iter+1) = lams(b, s, iter)*(1 - alpha) - beta * lams(exc(b), s, iter) + r(b, s, iter) + dw;
            plusNoise = thetaPairwiseDiffs;% + rand(N)*0.1;
            dtheta = rk4Step(@kuramoto, plusNoise, h, E{b}, lams(b, s, iter+1), omegas{b}, N(b));
            thetas{b}(:, s, iter+1) = thetas{b}(:, s, iter) + dtheta;

            % keep theta in [0, 2*pi)
            thetas{b}(:, s, iter+1) = wrapTo2Pi(thetas{b}(:, s, iter+1));

            % network mean phase vector
            z(b, s, iter+1) = sum(exp(1i*thetas{b}(:, s, iter+1))) / N(b);  
            r(b, s, iter+1) = abs(z(b, s, iter+1));
        end

        for b = 1:nNets
            if r(b, s, iter+1) > thres
                decision(b, s, :) = [b iter+1];
            end
        end

        if any(decision(:, s, :))
            return
        end
    end
end

function A_w = weightedA(net)
    % weighted adjacency table
    weightedEdges = table2array(net.Edges);
    % convert to matrix
    A_w = zeros(net.numnodes);
    for e = 1:size(weightedEdges, 1)
        A_w(weightedEdges(e, 1), weightedEdges(e, 2)) = weightedEdges(e, 3);
    end
end

function dTh = kuramoto(x, varargin)
    % Kuramoto oscillator differential equation 
    E_b = varargin{1};
    lam = varargin{2};
    omega = varargin{3};
    N_b = varargin{4};
    dTh = omega + (lam / N_b) * dot(E_b, sin(x))';
end

function res = rk4Step(dxdt, x, h, varargin)
    % 4-th order Runge-Kutta method on non-zero values in x
    f1 = dxdt(x, varargin{:});
    f2 = dxdt(x + 0.5*h*f1(:, 1), varargin{:});
    f3 = dxdt(x + 0.5*h*f2(:, 1), varargin{:});
    f4 = dxdt(x + h*f3(:, 1), varargin{:});        
    res = (h/6)*f1(:, 1) + 2*f2(:, 1) + 2*f3(:, 1) + f4(:, 1);
end
end