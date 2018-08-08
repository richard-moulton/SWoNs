function [theta, z] = kuramNetwork(net, lam, omega, theta0, steps, h)
% This script  runs a kuramoto oscillator through a pre-generated network
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
% BC/ML/SWoNS/2018
% Code adapted from appmath.wordpress.com, courtesy
% Jeongho Kim, Mathematical Sciences, Seoul National University.

% network parameters
N = net.numnodes;
E = weightedA(net);  % transposed because of `dot` in kuramoto

% pre-allocate arrays and set initial conditions
z = zeros(length(steps-1), 1);
theta = zeros(N, steps);
theta(:, 1) = theta0;

% partial function for evaluating the numerical integral
kuramotoPartial = @(thetaDiffs, varargin) kuramoto(thetaDiffs, omega, lam, N, varargin{:});

for iter = 1:steps-1
    % pairwise inter-node phase differences
    thetaPairwiseDiffs = theta(:, iter) - theta(:, iter)'; 
    
    % numerical integration step
    plusNoise = thetaPairwiseDiffs;% + rand(N)*0.1;
    dtheta = rk4Step(kuramotoPartial, plusNoise, h, ones(N));
    theta(:, iter+1) = theta(:, iter) + dtheta;
    
    % keep theta in [0, 2*pi)
    theta(:, iter+1) = wrapTo2Pi(theta(:, iter+1));
    
    % network mean phase vector
    z(iter+1) = sum(exp(1i*theta(:, iter+1))) / N;  
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

function dTh = kuramoto(x, omega, lam, N, varargin)
    % Kuramoto oscillator differential equation 
    E = varargin{1};
    dTh = omega + (lam / N) * dot(E, sin(x))';
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