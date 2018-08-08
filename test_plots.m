    r(iter+1) = abs (z);
    psi(iter+1) = angle(z);

figure;
for j = 1:length(h)
   subplot(ceil(length(h)/2),2,j)
   imagesc(log(1 + r(:,1:80,length(h))' - r(:,1:80,j)'))
   title(sprintf("h = %f", h(j)));
end

figure;
imagesc(r')
yticklabels(linspace(0,4,10))
xlabel('time')
ylabel('Coupling Constant')


stds = squeeze(std(theta2,[], 1))';
figure;
imagesc(stds);

figure;
imagesc(log(r' ./ stds(:, 1:end-1)))

%%
theta_i = theta{1, end};

figure;
for i = 1:size(theta_i, 2)
    
        x=cos(theta_i(:, i)');% - psi(i, :));
        y=sin(theta_i(:, i)');% - psi(i, :));
        
        s=linspace(0,2*pi,100);
        cx=cos(s);
        cy=sin(s);
        p = plot(x,y,'.',cx,cy);
        set(p,'MarkerSize',20);
        axis([-1 1 -1 1])
        axis square
        drawnow
        pause(0.1);
end
%%
figure;
r = abs(squeeze(z(1, :, :)));
imagesc(r)

%%
theta_p = theta{1, end};
z_p = squeeze(z(1, end, :));
r=abs(z_p);
psi=angle(z_p);

s=linspace(0,2*pi,100);
cx=cos(s);
cy=sin(s);

figure;
%p = plot(x,y,'.',cx,cy);
p=animatedline('Marker', '.', 'MarkerSize',20);
%p_z = plot([0 zx],[0 zy]); 
p_z=animatedline('LineWidth', 2);
%set(p_z, 'LineWidth', 2);
%set(p,'MarkerSize',20);
axis([-1 1 -1 1])
axis square

for i = 1:size(theta_p, 2)
        x=cos(theta_p(:, i)');% - psi(i, :));
        y=sin(theta_p(:, i)');% - psi(i, :));
        %set(p, 'XData', x);
        %set(p, 'YData', y);
        clearpoints(p);
        addpoints(p, x, y);
        
        
        zx = r(i) * cos(psi(i));
        zy = r(i) * sin(psi(i));
        %set(p_z, 'XData', zx);
        %set(p_z, 'YData', zy);  
        clearpoints(p_z);
        addpoints(p_z, zx, zy)
        drawnow 
        pause(0.1);
end

