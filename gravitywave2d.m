clc;clear all; close all

%% Set problem
xmax=500;  % x-direction water tank size
ymax=500;  % y-direction water tank size
dt=0.4;    % time step
dx=10;     % x-direction spatial step
dy=10;     % y-direction spatial step
g=9.8;     % gravity constant
% Create grid
xh=0:dx:xmax;
yh=0:dy:ymax;
[yh xh]=meshgrid(yh,xh);
u=zeros([size(xh) 2]); % initial u velocity
v=zeros([size(xh) 2]); % initial v velocity
eta=zeros([size(xh) 2]); % water surface height
h=ones([size(xh) 2])*20; % depth
% initial condition of water surface height
sigma=20;
eta(:,:,1)=1.5*exp(-(xh-xmax/2).^2/(2*sigma.^2)+-(yh-ymax/2).^2/(2*sigma.^2));
%% Solve 2d wave equation
for t=1:100
    
    % Solve momentum equation
    for i=2:size(xh,1)-1
        for j=2:size(yh,2)-1
            u(i,j,2)=u(i,j,1)-dt*g*(eta(i+1,j,t)-eta(i,j,t))/dx;
            v(i,j,2)=v(i,j,1)-dt*g*(eta(i,j+1,t)-eta(i,j,t))/dy;
        end
    end
        
    % Western boundary condition
    u(1,:)=0;
    v(1,:)=0;
    % Eastern boundary condition
    u(end-1,:)=0;
    v(end-1,:)=0;
    % Southern boundary condition
    u(:,1)=0;
    v(:,1)=0;
    % Northen boundary condition
    u(:,end-1)=0;
    v(:,end-1)=0;
    
    % Solve continuity equation
    for i=2:size(xh,1)-1
        for j=2:size(yh,2)-1
            up1=0.5*(u(i,j,2)+abs(u(i,j,2)));
            um1=0.5*(u(i,j,2)-abs(u(i,j,2)));
            vp1=0.5*(v(i,j,2)+abs(v(i,j,2)));
            vm1=0.5*(v(i,j,2)-abs(v(i,j,2)));
            up2=0.5*(u(i-1,j,2)+abs(u(i-1,j,2)));
            um2=0.5*(u(i-1,j,2)-abs(u(i-1,j,2)));
            vp2=0.5*(v(i,j-1,2)+abs(v(i,j-1,2)));
            vm2=0.5*(v(i,j-1,2)-abs(v(i,j-1,2)));
            
            eta(i,j,t+1)=eta(i,j,t)...
                -dt/dx*(up1*h(i,j)+um1*h(i+1,j)...
                                 -up2*h(i-1,j)-um2*h(i,j))...
                -dt/dy*(vp1*h(i,j)+vm1*h(i,j+1)...
                                 -vp2*h(i,j-1)-vm2*h(i,j));
        end
    end
    u(:,:,1)=u(:,:,2);
    v(:,:,1)=v(:,:,2);
    u(:,:,2)=0;
    v(:,:,2)=0;
end
figure
set(gcf,'color','w');
set(gcf,'Renderer','Zbuffer');
%% Visualization
for t=1:100
    surf(xh,yh,eta(:,:,t))
    zlim([-1 1])
    caxis([-1 1])
    drawnow
    title(['2d wave equation, t=' num2str((t-1)*dt,'%4.0f')]...
        ,'fontweight','bold');
    xlabel('x (m)','fontweight','bold')
    ylabel('y (m)','fontweight','bold')
    zlabel('\eta, sea surface height (m)','fontweight','bold')
    
    
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if t == 1;
        imwrite(imind,cm,'gravitywave2d.gif','gif','DelayTime',0.1, 'Loopcount',inf);
    else
        imwrite(imind,cm,'gravitywave2d.gif','gif','DelayTime',0.1,'WriteMode','append');
    end
end