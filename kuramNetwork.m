%% Run Kuramoto oscillator through network
clear; clc; clf; close all
% This script generates a network and runs a kuramoto oscillator
% conditioned on local connections
% Numerical integration through fourth-order Runge-Kutta Method
%% Generate network

N = 50;
K = 2;
q = 0.1;
h = 0.1;
Lam = 2;
displayFlag = false;
testNet = createNetwork (N, K, q,displayFlag);
Edges = testNet.Edges;
iter = 500;
numNodes = numel(unique(Edges.EndNodes(:,1)));
omega = randn(1,numNodes)*.1;  %initialize nodes with random intrinsic frequency
% omega(1:N/2) = omega(1:N/2) + 3;
theta = [2*pi*randn(1,numNodes)',zeros(numNodes,iter-1)];
f=figure;
reset(f); %Uncomment to save video
v = VideoWriter('Synchrony Unimodal.avi');
open(v)

for k = 1:iter-1
    
    thetaConnect = theta(:,k)*ones(1,N)-(ones(N,1)*theta(:,k)'); %Generates phase matrix
    
    for i = 1:numNodes
        
        tmp = Edges.EndNodes(find(Edges.EndNodes(:,1)==i),2);
        indEdgeConnect(i,:) = ismember(1:N,tmp); %generates connection matrix
        
    end
    
        thetaConnect(~indEdgeConnect) = 0; %no connection == 0
        
        f1 = kuramoto (thetaConnect,Lam,N,omega);
         f2=kuramoto(thetaConnect+0.5*h*f1(:,1),Lam,N,omega);
         f3=kuramoto(thetaConnect+0.5*h*f2(:,1),Lam,N,omega);           %4-th order Runge-Kutta method.
         f4=kuramoto(thetaConnect+h*f3(:,1),Lam,N,omega);        
        theta(:,k+1) = theta(:,k)+(h/6)*(f1(:,1))+2*f2(:,1)+2*f3(:,1)+f4(:,1);
        
        
        x=cos(theta(:,k));
        y=sin(theta(:,k));
        s=linspace(0,2*pi,100);
        cx=cos(s);
        cy=sin(s);
        p = plot(x,y,'.',cx,cy);
        set(p,'MarkerSize',20);
        axis([-1 1 -1 1])
        axis square
        A = getframe; %Uncomment to save video
        drawnow
        writeVideo(v,A)
    
end
% close(v)

function f=kuramoto(x,Lam,N,omega)
 
f=omega+(Lam/N)*sum(sin(x))';
 
end