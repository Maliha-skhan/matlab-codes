%%%%%%%%%%%%%%%%%%%%%% Start of the code %%%%%%%%%%%%%%%%%%%%%%
 clc
 clear all
 close all
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                                                                   %
 % Initial Condition                                                 %
 %                                                                   %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nx = 21; % number of nodes in x-direction
 ny = 21; % number of nodes in y-direction
 nt = 10; % number of time steps
 nit = 10; % number of artificial time steps
 % ...... for pressure
 vis = 0.1; % Viscosity
 rho = 1.0; % Density
 Lx = 2; % Length in x
 Ly = 2; % Length in y-direction
 
 dx = Lx/(nx-1); % grid spacing in x-direction
 dy = Ly/(ny-1); % grid spacing in y-direction
 dt = 0.01; % time step size
 x = 0:dx:Lx; % node x-ordinates
 y = 0:dy:Ly; % node y-ordinates
 u = zeros(ny,nx); % nodal velocity x-c
 v = zeros(ny,nx); % nodal velocity y-component
 p = zeros(ny,nx); % nodal pressure
 un = zeros(ny,nx); % time marched velocity x-dire
 vn = zeros(ny,nx); % time marched velocity x-direction
 pn = zeros(ny,nx); % temporary pressure for calculation

 b = zeros(ny,nx); % nodal source term value from pressure
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Boundary Condition                                                %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [y,x] = meshgrid(y,x);
 u(ny,:) = 1;

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                                                                   %
 % Calculation                                                       %
 %                                                                   %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

uT(:,:,1) = u(:,:); 
vT(:,:,1) = v(:,:);
pT(:,:,1) = p(:,:);

 for it = 1:nt+1 % loop over time
 for i = 2:nx-1 % this i,j loop is
 for j = 2:ny-1 % ... body force from pressure equation
     b(j,i) = rho*(((((u(j,i+1)-u(j,i-1))/(2*dx))+...
         ((v(j+1,i)-v(j-1,i))/(2*dy)))*(1/dt))+((u(j,i+1)...
         - u(j,i-1))/(2*dx)).^2+(2*((u(j,i+1)...
 -u(j,i-1))/(2*dy))*((v(j+1,i)-v(j-1,i))/(2*dx)))...
 +((v(j+1,i)-v(j-1,i))/(2*dy)).^2);
 end 
 end
 for iit = 1:nit+1
     pn=p;
     for i = 2:nx-1
for j = 2:ny-1
    p(j,i) = ((p(j,i+1)+p(j,i-1))*(dy^2)...
        +(p(j+1,i)+p(j-1,i))*(dx^2))/(2*(dx^2+dy^2))...
        +(b(j,i)*dx^2*dy^2)/(2*(dx^2+dy^2));
end
     end
     p(:,1) = p(:,2); p(:,nx) = p(:,nx-1);
     p(1,:) = p(2,:); p(ny,:) = p(ny-1,:);
 end
 
 un=u; vn = v; % assigning velocity values to nth-time values 
 % to calculate (n+1)th-time values
 for  i = 2:nx-1 % this i,j loop is for calculating velocity
 for j = 2:ny-1
     u(j,i) = un(j,i)-un(j,i)*(un(j,i)-un(j,i-1))*(dt/dx)...
-vn(j,i)*(un(j,i)-un(j-1,i))*(dt/dy)...
-(1/rho)*(p(j,i+1)-p(j,i-1))*(dt/(2*dx))...
+(vis/rho)*((un(j,i+1)-(2*un(j,i))+un(j,i-1))*...
(dt/dx^2)+(un(j+1,i)-(2*un(j,i))+un(j-1,i))*(dt/dy^2));

v(j,i) = vn(j,i)-un(j,i)*(vn(j,i)-vn(j,i-1))*(dt/dx)...
    -vn(j,i)*(vn(j,i)-vn(j-1,i))*(dt/dy)...
    -(1/rho)*(p(j+1,i)-p(j-1,i))*(dt/(2*dy))...
+(vis/rho)*((vn(j,i+1)-(2*vn(j,i))+vn(j,i-1))*...
(dt/dx^2)+(vn(j+1,i)-(2*vn(j,i))+vn(j-1,i))*(dt/dy^2));
 end
 end

 %%% Surface Plotting %%%%%%%%%%%%%%%%%%%%%% 
 figure(1)
 surf(x,y,u')
 title('Surface Plot: Velocity Component, u');
xlabel('x - ordinate');
ylabel('y - ordinate');
zlabel('Velocity component, u');

 figure(2)
 surf(x,y,v')
 title('Surface Plot: Velocity Component, v');
xlabel('x - ordinate');
ylabel('y - ordinate');
zlabel('Velocity component, v');

figure(3)
 surf(x,y,sqrt(u.^2+v.^2)')
 title('Surface Plot: Velocity field, U');
xlabel('x - ordinate');
ylabel('y - ordinate');
zlabel('Velocity component, v');

figure(4)
 surf(x,y,p')
 title('Surface Plot: Pressure, P');
xlabel('x - ordinate');
ylabel('y - ordinate');
zlabel('Pressure, P');

%%%%%%%%%%%%%%%%%%%%%% Contour Plotting %%%%%%%%%%%%%%%%%%%%%%
figure(5)
contourf(x,y,u',10)
 title('Contour: Velocity component, u');
xlabel('x - ordinate');
ylabel('y - ordinate');

figure(6)
contourf(x,y,v',10)
 title('Contour: Velocity component, v');
xlabel('x - ordinate');
ylabel('y - ordinate');

figure(7)
contourf(x,y,sqrt(u.^2+v.^2)',10)
 title('Contour Plot: Velocity field, U');
xlabel('x - ordinate');
ylabel('y - ordinate');


figure(8)
contourf(x,y,p',10)
 title('Contour Plot: Pressure, P');
xlabel('x - ordinate');
ylabel('y - ordinate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(:,1) = 0; u(:,nx) = 0; u(1,:) = 0; u(ny,:) = 1; 
v(:,1) = 0; v(:,nx) = 0; v(1,:) = 0; v(ny,:) = 0;
 end
 %%%%%%%%%%%%%%%%%%%%%%% End of the code %%%%%%%%%%%%%%%%%%%%%%% 
 