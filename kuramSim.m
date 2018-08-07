clc, clf, clear, close all
 
% Numerical simulation for the Kuramoto model:
  
N=50;   % number of node.
K=1;    % Coupling strength.
h=0.1;
iter=500;
t=0:h:h*iter;
theta=zeros(N,iter);
theta(:,1)=2*pi*rand(N,1);
figure;
subplot(2,1,1)
Omega(1:N/2) = .1;
Omega(N/2:N) = .9;
Omega = Omega';
hist(Omega,50);
title('Distribution of Frequencies');

for j=1:iter
 
k1=kuramoto(theta(:,j),K,N,Omega);
theta(:, j+1)=theta(:,j)+(h/6)*(k1+2);%*k2+2*k3+k4);
x=cos(theta(:,j));
y=sin(theta(:,j));
s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);
subplot(2,1,2)
plot(x,y,'o',cx,cy)
axis([-1 1 -1 1])
axis square
 
drawnow
 
end
title('Phase')