clear all; close all; clc;
%% Constants
G = 0.01; % gravitational constant
F = 0; % coriolis frequency
BETA = 2E-11; % beta-plane parameter
T1 = 200 * 24 * 3600; % end time (200 days)
D = 400; % water depth
X = 2000e3; % domain size
Y = 4000e3;
KAPPA = 1/(100 * 24 * 3600); % bottom friction parameter (1/8.5 days)
LAMBDA = 1/(200 * 24 * 3600); % damping coefficient (1/200 days)
WIDTH = 400e3; % source width
Q = 10e6; % source strength (10 Sv)
DX = 40e3; % target grid spacing
DY = 40e3;
DT_FACTOR = 10; % you may have to increase this to ensure stability

%% Construct staggered grid
x = linspace(0,X,2*fix(X/DX));
y = linspace(-Y/2,Y/2,2*fix(Y/DY));

[X_H, Y_H] = meshgrid(x(1:2:end),y(1:2:end));
[X_U, Y_U] = meshgrid(x(2:2:end),y(1:2:end));
[X_V, Y_V] = meshgrid(x(1:2:end),y(2:2:end));

DX = 2*(x(2) - x(1))
DY = 2*(y(2) - y(1))
DT = min(DX,DY) / sqrt(G*D) / DT_FACTOR % calculate dt from CFL condition

%% Initial conditions
% rest
h0 = ones(size(X_H)) * D;
u0 = zeros(size(X_U));
v0 = zeros(size(X_V));

figure(1)
clf();

subplot(3,1,1);
surf(X_H,Y_H,h0);
xlabel('X');
ylabel('Y');
zlabel('Layer height');
title('Initial conditions');

subplot(3,1,2);
surf(X_U,Y_U,u0);
xlabel('X');
ylabel('Y');
zlabel('U');

subplot(3,1,3);
surf(X_V,Y_V,v0);
xlabel('X');
ylabel('Y');
zlabel('V');

%% Prepare source term
% require that the integrated source is equal to Q
water_source = exp(-(X_H.^2 + (Y_H-Y/4).^2)/WIDTH^2);
source_normalization = sum(reshape(water_source,1,[])*DX*DY);
water_source = Q / source_normalization * water_source;

figure(5);
clf();
contourf(X_H,Y_H,water_source);
xlabel('X');
ylabel('Y');
colorbar();
title('Source strength (m^3/s)');

%% Prepare output
figure(2);
clf();
%plot_h = surf(X_H,Y_H,h0);
[~, plot_h] = contourf(X_U,Y_U,h0,390:5:475);
xlabel('X');
ylabel('Y');
zlabel('Layer height');
title_h = title('t = 0');
c = colorbar();
caxis([400 450])

figure(3);
clf();
colormap redbluecmap
[~, plot_u] = contourf(X_U,Y_U,u0);
xlabel('X');
ylabel('Y');
title_u = title('t = 0');
c = colorbar();
c.Label.String = 'U';

figure(4);
clf();
colormap redbluecmap
[~, plot_v] = contourf(X_V,Y_V,v0);
xlabel('X');
ylabel('Y');
title_v = title('t = 0');
c = colorbar();
c.Label.String = 'V';

figure(2);

%% Main solution loop
t = 0;
h = h0;
u = u0;
v = v0;

while t < T1
    %% update solution
    
    % To avoid loops and explicit treatment of the boundaries, we introduce
    % ghost grid cells that represent the boundary conditions for u and h.
    
    % add ghost cells to make h large enough to shift it by 1 and -1
    % note that h_pad(2:end-1,2:end-1) = h
    h_pad = [nan, h(1,:), nan ; h(:,1), h, h(:,end) ; nan, h(end,:), nan];
    
    % update solution for u and v
    % h padding leads to u_new(:,end) == v_new(end,:) == 0 for our
    % staggered grid
    
    % calculate prediction
    u_pred = (1-DT*KAPPA) * u - G*DT/DX*(h_pad(2:end-1,3:end) - h);
    v_pred = (1-DT*KAPPA) * v - G*DT/DY*(h_pad(3:end,2:end-1) - h);
    
    % update with coriolis force
    alpha = DT * (F + BETA * Y_H);
    beta = 0.25*alpha.^2;
    u_new = (u_pred - beta .* u + alpha .* v)./(1+beta);
    v_new = (v_pred - beta .* v - alpha .* u)./(1+beta);
    
    % re-enforce BC
    u_new(:,end) = 0;
    v_new(end,:) = 0;
    
    % add ghost cells to u_new and v_new (no-slip BC)
    [m, n] = size(u_new);
    u_new_pad = [zeros(1,n+2) ; zeros(m,1), u_new, zeros(m,1) ; zeros(1,n+2)];
    u_plus = .5*(u_new_pad + abs(u_new_pad));
    u_minus = .5*(u_new_pad - abs(u_new_pad));
    
    [m, n] = size(v_new);
    v_new_pad = [zeros(1,n+2) ; zeros(m,1), v_new, zeros(m,1) ; zeros(1,n+2)];
    v_plus = .5*(v_new_pad + abs(v_new_pad));
    v_minus = .5*(v_new_pad - abs(v_new_pad));
    
    % calculate upwinding values
    uh_e = u_plus(2:end-1,2:end-1) .* h + u_minus(2:end-1,2:end-1) .* h_pad(2:end-1,3:end);
    vh_n = v_plus(2:end-1,2:end-1) .* h + v_minus(2:end-1,2:end-1) .* h_pad(3:end,2:end-1);
    uh_w = u_plus(2:end-1,1:end-2) .* h_pad(2:end-1,1:end-2) + u_minus(2:end-1,1:end-2) .* h;
    vh_s = v_plus(1:end-2,2:end-1) .* h_pad(1:end-2,2:end-1) + v_minus(1:end-2,2:end-1) .* h;
    
    % update solution for h
    h_new = h - DT/DX*(uh_e - uh_w) - DT/DY*(vh_n - vh_s) + DT*(water_source - LAMBDA*(h-D));

    %% update solutions
    h = h_new;
    u = u_new;
    v = v_new;
    
    %% step time
    t = t + DT;
    
    %% output
    plot_h.ZData = h;
    title_h.String = ['t = ' num2str(t)];
    plot_u.ZData = u;
    title_u.String = ['t = ' num2str(t)];
    plot_v.ZData = v;
    title_v.String = ['t = ' num2str(t)];
    drawnow limitrate;
    
end
% Conclusion: KAPPA=1/(100 days)
% Strong difussion coefficient 
% more upwelling everywhere 