clear all; close all; clc;
H_ice=3027.3; %in m
H=H_ice-25;
d=[102.5 202.5 302.5 402.5 502.5 602.5 702.5 767.5 872.5 1022.5 1221.25 1391.25];
%lambda=[0.256 0.213 0.204 0.204 0.182 0.192 0.169 (5/26) (5/34) (5/30) (2.5/(19+1.5))];
lambda2=[5/21.5 5/26 5/24.5 5/24.5 5/28.5 5/26 5/31 (5/27.5) (5/34) (5/31) (2.5/(22)) (2.5/25.5)];
%lambda2=[5/22 5/26.5 5/25.5 5/25.5 5/28.5 5/27.5 5/32.5 (5/29.5) (5/35.5) (5/32) (2.5/(22.5)) (2.5/26)];
z=H_ice-d;

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( z, lambda2 );

% Set up fittype and options.
ft = fittype( 'p1*x+p2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.672400066590172 0.930103078748882];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'Data', 'Linear fit', 'Location', 'NorthEast' );
% Label axes
ylabel('\lambda [m]')
xlabel('z [m]')
grid on


% Bestem Dansgaard-Johnsen model-parametrene
p1=8.786e-05; p2=-0.04064;
%p1 =9.759e+05; p2 =674.5;

b=((p1*H_ice)+p2); %akkumulationen i noterne er det "c"
h0=H_ice % Is tykkelse
h1=2*p2; % højde hvor knækket befinder sig. Her er u ikke længere konstant. I noterne "h'"   

% time function
a=( (2*h0)-h1 )/b; %hældning
t1 = (1/2) * a * log( ( (2*h0) - h1 )/h1 )
z1=(h1/50):(h1/50):h1;
t = ( a*( (h1./z1) - 1 ) ) + t1;
%t=a*((h1./z)-1);
figure(2); hold on; grid on; box on;
plot(z1,t,'bd'); ylabel('time'); xlabel('depth (m)')

% Fit: 'untitled fit 1'.
[xData2, yData2] = prepareCurveData( z1, t );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult2, gof2] = fit( xData2, yData2, ft, 'Normalize', 'on' );


% vulkan udbrud 1
%z=40.3
t_vulkan1 =( a*( (h1/40.3) - 1 ) )+ t1
t_vulkan2 =( a*( (h1/58.5) - 1 ) )+ t1
t_vulkan3 =( a*( (h1/67.1) - 1 ) )+ t1
t_vulkan4 =( a*( (h1/187.3) - 1 ) )+ t1
t_vulkan5 =( a*( (h1/256.3) - 1 ) )+ t1 
t_vulkan6 =( a*( (h1/736.70 ) - 1 ) )+ t1 