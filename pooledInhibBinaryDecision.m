function [thetas, r, lams, decision] = pooledInhibBinaryDecision(nets, thres, alpha, beta, eps, omegas, thetas0, lams0, maxsteps, h)
% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.

% network parameters
for b = 1:size(nets, 2)
    N(b) = nets{b}.numnodes;
    E{b} = weightedA(nets{b});
    
    % pre-allocate arrays and set initial conditions
    thetas{b} = zeros(N(b), maxsteps);
    thetas{b}(:, 1) = thetas0{b};
    lams(b, 1) = lams0(b);
    z(b, 1) = sum(exp(1i*thetas{b}(:, 1))) / N(b); 
    r(b, 1) = abs(z(b, 1)); 
end

% initialize decision matrix
decision = zeros(size(nets, 2), 2);

exc = [2 1];

for iter = 1:maxsteps-1
    for b = 1:size(nets, 2)
        % pairwise inter-node phase differences
        thetaPairwiseDiffs = thetas{b}(:, iter) - thetas{b}(:, iter)'; 

        % numerical integration step
        dw = rand * eps;
        lams(b, iter+1) = lams(b, iter)*(1 - alpha) - beta * lams(exc(b), iter) + r(b, iter) + dw;
        %* log(r(b, iter) / sum(r(:, iter)));  
        plusNoise = thetaPairwiseDiffs;% + rand(N)*0.1;
        dtheta = rk4Step(@kuramoto, plusNoise, h, E{b}, lams(b, iter+1), omegas{b}, N(b));
        thetas{b}(:, iter+1) = thetas{b}(:, iter) + dtheta;

        % keep theta in [0, 2*pi)
        thetas{b}(:, iter+1) = wrapTo2Pi(thetas{b}(:, iter+1));

        % network mean phase vector
        z(b, iter+1) = sum(exp(1i*thetas{b}(:, iter+1))) / N(b);  
        r(b, iter+1) = abs(z(b, iter+1));
    end
    
    for b = 1:size(nets, 2)
        if r(b, iter+1) > thres
            decision(b, :) = [1 iter+1];
        end
    end
    
    if any(decision)
        return
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