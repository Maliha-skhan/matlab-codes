clear all; close all; clc; 
% Set problem
xmax=2000*10^3;  % x-direction
%ymax=2000*10^3;  % y-direction
nmax=300; % max time steps
dt=720;    % time step
dx=20*10^3;     % x-direction spatial step
dy=dx;     % y-direction spatial step
ymax=1000*10^3;
%ymax = linspace(-L,L,101)'; % Coordinate vector
%g=9.8;     % gravity constant
g=.01; 

beta=10*2*10^(-11); %Rossby beta 
f0=10^(-4);
%f0=0;
a=dt*f0; %alpha
b=(1/4)*a.^2; %beta

% Create grid
xh=0:dx:xmax; xplot = xh;
yh=-ymax:dy:ymax; y=yh;
[yh xh]=meshgrid(yh,xh);
%f0=f*[];
u=zeros([size(xh) 2]); % initial u velocity
v=zeros([size(yh) 2]); % initial v velocity
eta=zeros([size(xh) 2]); % water surface height
D=500; %depth
%lambda=180*10^3;%;72000; 
lambda=180*10^3;%;72000; 
       k=2*pi/lambda;
       % Eps=1-( sin(0.5*k*dx) /(0.5*k*dx) );  
   Eps=zeros([size(xh)]);
       %Eps=1-( sin(0.5*k*dy) /(0.5*k*dy) );   
       Eps(1:20,:)=0.5;
       Eps(81:101,:)=0.5;
              Eps(:,1:20)=0.5;
       Eps(:,81:101)=0.5;
       %f=zeros([101 101])     
       %f(:,1:101)=f0+(beta.*y(1:101))';
      %& f(:,:)=f(1:101)';
% initial condition of water surface height
sigma=130*10^3;


eta(:,:,1)=0.8*exp(-(xh-xmax/2).^2/(sigma.^2)+-(yh).^2/(sigma.^2));
for j=1:101
                %f(1:101,j)=f0+(beta*y(j));
                %f(1:101,51)=(beta*y(51));
                f(1:101,j)=beta*y(j);
end
a=dt*f(:,:); %alpha
b=(1/4)*a.^2; %beta
%%
% Solve 2d wave equation
 %f(1:101,1)=f0+(beta*y(1));
for n=1:nmax
    h(:,:,n)=eta(:,:,n)+D;
    % Solve momentum equation
    for i=2:size(xh,1)-1
        for j=2:size(yh,2)-1
%            f(1:101,j-1)=f0+(beta*y(j));
            us(i,j,n+1)=u(i,j,n)-dt/dx*g*(eta(i+1,j,n)-eta(i,j,n));
            vs(i,j,n+1)=v(i,j,n)-dt/dy*g*(eta(i,j+1,n)-eta(i,j,n));
        %%% Now we add the coriolis parameter
                    u(i,j,n+1)=( us(i,j,n+1) - b(i,j)*u(i,j,n)+a(i,j)*v(i,j,n) )/(1+b(i,j));
            v(i,j,n+1)=( vs(i,j,n+1)-b(i,j)*v(i,j,n)-a(i,j)*u(i,j,n) )/(1+b(i,j));
                        
        end
    end
          
    % Eastern boundary condition
    u(end,:,:)=0;
    us(end,:,:)=0;    
    
    % Northen boundary condition
    v(:,end,:)=0;
    vs(:,end,:)=0;
    
    % Solve continuity equation
    for i=2:size(xh,1)-1
        for j=2:size(yh,2)-1
            up1=0.5*(u(i,j,n+1)+abs(u(i,j,n+1)));
            um1=0.5*(u(i,j,n+1)-abs(u(i,j,n+1)));
            
            vp1=0.5*(v(i,j,n+1)+abs(v(i,j,n+1)));
            vm1=0.5*(v(i,j,n+1)-abs(v(i,j,n+1)));
            
            up2=0.5*(u(i-1,j,n+1)+abs(u(i-1,j,n+1)));
            um2=0.5*(u(i-1,j,n+1)-abs(u(i-1,j,n+1)));
            
            vp2=0.5*(v(i,j-1,n+1)+abs(v(i,j-1,n+1)));
            vm2=0.5*(v(i,j-1,n+1)-abs(v(i,j-1,n+1)));


            eta(i,j,n+1)=eta(i,j,n)-(dt/dx*(up1*h(i,j)+um1*h(i+1,j)-up2*h(i-1,j)-um2*h(i,j)))-(dt/dy*(vp1*h(i,j)+vm1*h(i,j+1)-vp2*h(i,j-1)-vm2*h(i,j)));
       %Shapiro filter
%       k=2*pi/180000;
%       Eps=1-(sin(k*dx/2)/(k*dx/2));    
       %eta(i,j,n+1)=(1-Eps)*eta(i,j,n+1) + 0.25*Eps*(eta(i+1,j+1,n+1)+eta(i-1,j-1,n+1));
       %eta(i,j,n+1)=(1-Eps(i,j))*eta(i,j,n+1) + 0.25*Eps(i,j)*(eta(i+1,j,n+1)+eta(i-1,j,n+1));
       %eta(i,j)=(1-Eps)*eta(i,j) + 0.25*Eps*(eta(i+1,j)+eta(i-1,j));   

end
    end     

end

%%%Prepare videoobject
%writerObj = VideoWriter('GravityWaveCoriolis2D.avi');
%writerObj.FrameRate = 14; % How many frames per second.
%open(writerObj);

%%figur
%figure(1); %hold on
set(gcf,'color','w');
set(gcf,'Renderer','Zbuffer');

% Visualization
for n=1:nmax
    surf(xh,yh,eta(:,:,n),'FaceAlpha',0.9);colorbar; colormap jet;
    zlim([-0.3 0.8])
    caxis([-0.5 0.5])
    title(['Equatorial Waves, [hr] t_{hr} = ' num2str((1/3600)*(n)*dt) ', [days] t_{day} = ' num2str((1/(24*3600))*(n)*dt) ] ,'fontweight','bold');
   xlabel('E-W (m)','fontweight','bold')
    ylabel('N-S (m)','fontweight','bold')
    zlabel('\eta, sea surface height (m)','fontweight','bold')
     drawnow

%%%  Save videoframes
  %frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
  %writeVideo(writerObj, frame);    
 
end
%%
%LR=250*10^3; %
LR=(sqrt(g*D))/f0;
c_theory=-(beta*(LR^2))/3; % equatorial phase speed for the Rossby Wave
%%
%%%Close videoobject
%close(writerObj); %closing and saving movie


%%
figure(11);
ts4=eta(:,48,round(150));  ts4=ts4(:,:)' ;
time=dt*(1:nmax+1);
plot(xplot,ts4); title(['Positions series for point for fixed y=' num2str((0.34*ymax/dy)) ,'and N=' num2str(nmax/2)]); xlabel('x'); ylabel('\eta : sea surface height'); grid on; box on;

figure(10);
ts4=eta(:,38,446);  ts4=ts4(:,:)' ;
time=dt*(1:nmax+1);
plot(xplot,ts4); title(['Positions series for point for fixed y = ' num2str((38*ymax/dy)) ,'m and t = ' num2str(6.2424e+04),'s']); xlabel('x'); ylabel('\eta : sea surface height'); grid on; box on;



figure(12);
ts4=eta(:,(0.38*ymax/dy),round(nmax*0.3));  ts4=ts4(:,:) ;
time=dt*(1:nmax+1);
plot(xplot,ts4); title(['Positions series for point for fixed y=' num2str((0.38*ymax/dy)) ,'and N=' num2str(nmax*0.3)]); xlabel('x'); ylabel('\eta : sea surface height'); grid on; box on;

figure(70);
ts4=eta(:,(25),round(nmax*0.75));  ts4=ts4(:,:) ;
time=dt*(1:nmax+1);
plot(xplot,ts4); title(['Positions series for point for fixed y=' num2str((25)) ,'and t=' num2str(nmax*dt*0.75),'s']); xlabel('x'); ylabel('\eta : sea surface height'); grid on; box on;


figure(13);
ts=eta((round(xmax/dx*0.36)),(round(ymax/dy*0.38)),:);  ts=ts(:,1:nmax+1)' ;
time=(dt*(1:nmax+1))/3600;
plot(time,ts); title(['Time series for point (x,y)= (' num2str(xmax/dx*0.36),',' num2str(xmax/dx*0.38),')']);
xlabel('time [hr]'); ylabel('\eta : sea surface height'); grid on; box on;


ts2=eta(:,(ymax/dy*0.33),:); ts2=ts2(:,:);
figure(14); contour(ts2) ; xlabel('N`th step in time') ; ylabel('i`th step in x direction');  
colorbar;  title(['\eta contour plot, y=' num2str(ymax/dx*0.33)]); grid on; box on


%ts2=eta(:,(ymax/dy*0.5),:); ts2=ts2(:,:);
%figure(17); contour(ts2) ; xlabel('N`th step in time') ; ylabel('i`th step in x direction');  
%colorbar;  title(['\eta contour plot, y=' num2str(ymax/dx*0.5)]); grid on; box on


figure(18);
ts=eta((round(xmax/dx*0.36)),(round(ymax/dy*0.78)),:);  ts=ts(:,1:nmax+1)' ;
time=(dt*(1:nmax+1))/3600;
plot(time,ts); title(['Time series for point (x,y)= (' num2str(xmax/dx*0.36),',' num2str(xmax/dx*0.78),')']);
xlabel('time [hr]'); ylabel('\eta : sea surface height'); grid on; box on;



%% Phase speed determination 
%lambda=;
%T=;
%c_p=lambda/T;
%%
figure(22);
    h=surf(xh,yh,eta(:,:,nmax/2),'FaceAlpha',0.7);colorbar; c = jet(20); colormap(c) %colormap winter;
    %set(h, 'edgecolor','none')
    zlim([-0.02 .03])
    caxis([-.015 .015])
    title(['Gravity waves in 2D, [hr] t = ' num2str((1/3600)*(nmax/2)*dt) ] ,'fontweight','bold');
   xlabel('E-W (m)','fontweight','bold')
    ylabel('N-S (m)','fontweight','bold')
    zlabel('\eta, sea surface height (m)','fontweight','bold')
     drawnow
     
     %%
figure(23);
    h=surf(xh,yh,eta(:,:,1),'FaceAlpha',0.9);colorbar; c = jet(20); colormap(c) %colormap winter;
    set(h, 'edgecolor','none'); grid on
    zlim([0 1])
    caxis([-.25 .5])
    title(['Equatorial waves at initial state, [hr] t = ' num2str((1/3600)*(1)*dt) ] ,'fontweight','bold');
   xlabel('E-W (m)','fontweight','bold')
    ylabel('N-S (m)','fontweight','bold')
    zlabel('\eta, sea surface height (m)','fontweight','bold')
     drawnow
          %%
figure(44);
    h=surf(xh,yh,eta(:,:,round((4/5*nmax))),'FaceAlpha',0.9); colorbar; c = jet(20); colormap(c) %colormap winter;
    set(h, 'edgecolor','none'); grid on
    zlim([-0.3 1])
    caxis([-.25 .5])
    title(['Equatorial waves, [hr] t = ' num2str((1/3600)*((4/5*nmax))*dt) ] ,'fontweight','bold');
   xlabel('E-W (m)','fontweight','bold')
    ylabel('N-S (m)','fontweight','bold')
    zlabel('\eta, sea surface height (m)','fontweight','bold')
     drawnow
     %%
     figure(24); eta11=eta(:,50,:); eta11=eta11(:,:);
     contourf(xplot,(round(time*60*60)),eta11'); ylabel('time [sec]'); xlabel('x position [m]');
     colorbar; caxis([-.150 .150]); %axis([min(time) max(time) min(xplot) max(xplot)])
     title('Contour plot for fixed y value at the middle of the grids')
          %%
     figure(25); eta11=eta(:,70,:); eta11=eta11(:,:);
     contourf(xplot,(round(time*60*60)),eta11'); ylabel('time [sec]'); xlabel('x position [m]');
     colorbar; caxis([-.150 .150]); %axis([min(time) max(time) min(xplot) max(xplot)])
     title('Contour plot for fixed y value at y=70 of the grids')
    
     %%
c_kelvin=((1.76-1.18)*10^6)/((3.5114*10^5)-(7.848*10^4)) % equatorial kelvin wave phase speed