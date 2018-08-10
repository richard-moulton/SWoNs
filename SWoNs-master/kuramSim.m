% clc, clf, clear, close all
 
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
Omega= randn(1,50);
Omega = Omega';
hist(Omega,50);
title('Distribution of Frequencies');
Lam = linspace(0,4,50);

for k = 1:length (Lam)
    for j=1:iter
 
        k1=kuramoto(theta(:,j),Lam(k),N,Omega);
        theta(:, j+1)=theta(:,j)+(h/6)*(k1+2);%*k2+2*k3+k4);
        indOver = find (theta(:,k+1) > 2*pi);
        indUnder = find (theta(:,k+1) < 0);

        theta(indOver,k+1) = theta(indOver,k+1) - 2*pi;
        theta(indUnder,k+1) = theta(indUnder,k+1) + 2*pi;
        
        z = sum(exp(1i*sin(theta(:,j+1))))/N;
        r(k,j+1) = abs (z);
        psi(k,j+1) = angle(z);

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
end
figure; imagesc(r');

col = brewermap (50,'*Rd');
figure; 
for i = 1:50
    hold on
    plot(r(i,:),'Color',col(i,:))
end