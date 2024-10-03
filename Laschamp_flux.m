clear all;close all;clc;
data=importdata('Data_Laschamp_CP.xlsx');
window_size=10;
depth_ngrip=data.data.NGRIP(:,1);
Be10_ngrip=data.data.NGRIP(:,2);

depth_NGRIP = [2106.01 2109.62  2111.58  2113.22  2115.41  2118.62  2129.82  2132.64];
GICC05_age = [40563 40794 40912 41002 41109 41249 41858 42067];

% figure(1); plot(depth_ngrip,Be10_ngrip,'b'); title('NGRIP'); hold on; xlabel('depth [m]'); 
% ylabel('^{10}Be concentration [10^{3} atoms/g]'); grid on; box on;
% hold off

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth_NGRIP, GICC05_age );

% Set up fittype and options.
ft = fittype( 'a*x+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.504851666205104 0.442315182701475];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'untitled fit 1' );
figure(11);
h = plot( fitresult, xData, yData );
legend( h, 'GICC05_age vs. depth_NGRIP', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth_NGRIP
ylabel GICC05_age
grid on; close figure 11;

figure(2);%hold on;
%subplot(2,1,1);
fit1=fitresult((depth_ngrip));
plot(fit1,Be10_ngrip,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be concentration [10^{3} atoms/g]');
hold on; %window_size = 12;
Be10_ngrip(isnan(Be10_ngrip))=0;
modif = tsmovavg(Be10_ngrip,'m',window_size,1);  
plot(fitresult((depth_ngrip)),modif,'k', 'Linewidth',1); hold on; 
grid on; box on;
   %title('NGRIP'); %
   title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off
%%
depth_vostok=data.data.Vostok(:,1);
Be10_vostok=data.data.Vostok(:,2);
%
depth_VOSTOK_5G = [598.16  600.75 603.99 614.86 617.40];
GICC05_age = [40563 40794 41002 41858 42067];
%
[xData, yData] = prepareCurveData( depth_VOSTOK_5G, GICC05_age );
% Set up fittype and options.
ft = fittype( 'm*x+n', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.957506835434298 0.964888535199277];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% Plot fit with data.
figure(12);
h = plot( fitresult, xData, yData );
legend( h, 'GICC05_age vs. depth_VOSTOK_5G', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth_VOSTOK_5G
ylabel GICC05_age
grid on; close figure 12;

figure(4);
%figure(2);hold on;
%subplot(2,1,2);
fit2=fitresult(depth_vostok);
plot(fit2,Be10_vostok,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be concentration [10^{3} atoms/g]');
hold on; %window_size = 12;
Be10_vostok(isnan(Be10_vostok))=0;
modif2 = tsmovavg(Be10_vostok,'m',window_size,1);  
plot(fitresult(depth_vostok),modif2,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('Vostok,GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off
%%
window_size=5;
depth_edc=data.data.EDC(:,1);
Be10_edc=data.data.EDC(:,2);
%
depth_EDC = [731.65 734.55   736.97  738.19  739.80  746.68  749.17];
GICC05_age = [40563 40794 41002 41109 41249 41858 42067];

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth_EDC, GICC05_age );

% Set up fittype and options.
ft = fittype( 'd*x+f', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.913375856139019 0.0975404049994095];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(13);
h = plot( fitresult, xData, yData );
legend( h, 'GICC05_age vs. depth_EDC', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth_EDC
ylabel GICC05_age
grid on; close figure 13

figure(5);
%figure(2);hold on;
%subplot(2,1,2);
fit3=fitresult(depth_edc);
plot(fit3,Be10_edc,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be concentration [10^{3} atoms/g]');
hold on; %window_size = 12;
Be10_edc(isnan(Be10_edc))=0;
modif3 = tsmovavg(Be10_edc,'m',window_size,1);  
plot(fitresult(depth_edc),modif3,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('EDC, GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off


%%
depth_edml=data.data.EDML(:,1);
Be10_edml=data.data.EDML(:,2);
%
depth_EDML = [1362.28 1366.56  1369.54  1370.57 1372.73   1375.15  1386.61 1390.49];
GICC05_age = [40563 40794 40912 41002 41109 41249 41858 42067];
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth_EDML, GICC05_age );

% Set up fittype and options.
ft = fittype( 'h*x+s', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.957166948242946 0.970592781760616];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure(14);
h = plot( fitresult, xData, yData );
legend( h, 'GICC05_age vs. depth_EDML', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth_EDML
ylabel GICC05_age
grid on
close figure 14

figure(6);
%figure(2);hold on;
%subplot(2,1,2);
fit4=fitresult(depth_edml);
plot(fit4,Be10_edml,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be concentration [10^{3} atoms/g]');
hold on; %window_size = 12;
Be10_edml(isnan(Be10_edml))=0;
modif4 = tsmovavg(Be10_edml,'m',window_size,1);  
plot(fitresult(depth_edml),modif4,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('EDML, GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off


%%
% %%
% figure(7); hold on; 
% subplot(4,1,1); hold on
% plot(fit1,Be10_ngrip,'b'); hold on; %xlabel('year b2k'); 
% title('^{10}Be conc., GICC05 timescale')
% %ylabel('^{10}Be conc. [10^{3} atoms/g]');
% hold on; window_size = 5;
% plot(fit1,modif,'k', 'Linewidth',1); hold on; 
% grid on; box on; axis([min(fit1) max(fit2) min(Be10_ngrip) max(Be10_ngrip)])
%    %title('NGRIP'); %
%    legend('NGRIP')
% %   legend('data','moving average')
%       h=gca; 
%    get(h,'FontSize') 
%    set(h,'FontSize',12);
% hold off
% 
% subplot(4,1,2); hold on
% plot(fit2,Be10_vostok,'b'); hold on; %xlabel('year b2k'); 
% %ylabel('^{10}Be conc. [10^{3} atoms/g]');
% hold on; window_size = 5;
% plot(fit2,modif2,'k', 'Linewidth',1); hold on; 
% grid on; box on; axis([min(fit1) max(fit2) min(Be10_vostok) max(Be10_vostok)])
%    %title('NGRIP'); %
%    legend('Vostok')
% %   legend('data','moving average')
%       h=gca; 
%    get(h,'FontSize') 
%    set(h,'FontSize',12);
% hold off
% 
% subplot(4,1,3); hold on
% plot(fit3,Be10_edc,'b'); hold on; %xlabel('year b2k'); 
% ylabel('^{10}Be conc. [10^{3} atoms/g]');
% hold on; window_size = 5;
% plot(fit3,modif3,'k', 'Linewidth',1); hold on; 
% grid on; box on; axis([min(fit1) max(fit2) min(Be10_edc) max(Be10_edc)])
%    %title('NGRIP'); %
%    legend('EDC')
% %   legend('data','moving average')
%       h=gca; 
%    get(h,'FontSize') 
%    set(h,'FontSize',12);
% hold off
% 
% subplot(4,1,4); hold on
% plot(fit4,Be10_edml,'b'); hold on; xlabel('year b2k'); 
% %ylabel('^{10}Be conc. [10^{3} atoms/g]');
% hold on; window_size = 5;
% plot(fit4,modif4,'k', 'Linewidth',1); hold on; 
% grid on; box on; axis([min(fit1) max(fit2) min(Be10_edml) max(Be10_edml)])
%    %title('NGRIP'); %
%    legend('EDML')
% %   legend('data','moving average')
%       h=gca; 
%    get(h,'FontSize') 
%    set(h,'FontSize',12);
%    %figure(7);hold on; ylabel('^{10}Be conc. [10^{3} atoms/g]'); hold on;
% hold off

%% See: ftp://ftp.ncdc.noaa.gov/pub/data/paleo/icecore/greenland/summit/ngrip/ngrip-10be.txt
%Density of ice: 917 kg/m^3

an=31556926; % 1 year = 31 556 926 seconds

d=917*10^3;
data2=importdata('AICC2012_official.xls');
acc_ngrip=data2.data.NGRIP(:,6); acc_ngrip=acc_ngrip.*(d/an);
acc_ngrip=acc_ngrip(2095:2133);
yr_ngrip=data2.data.NGRIP(:,2)-1950; yr_ngrip=yr_ngrip(2095:2133);
depth2NGRIP=data2.data.NGRIP(:,1); depth2NGRIP=depth2NGRIP(2095:2133);
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth2NGRIP, acc_ngrip );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure(19);
h = plot( fitresult, xData, yData );
legend( h, 'acc_ngrip vs. depth2NGRIP', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth2NGRIP
ylabel acc_ngrip
grid on; close figure 19;
ACC_NGRIP=fitresult(depth_ngrip); flux_ngrip=Be10_ngrip.*ACC_NGRIP*10^3;


%%
acc_edc=data2.data.EDC(:,6); acc_edc=acc_edc.*(d/an); acc_edc=acc_edc(424:478);
yr_edc=data2.data.EDC(:,2)-1950; yr_edc=yr_edc(424:478);
depth2EDC=data2.data.EDC(:,1); depth2EDC=depth2EDC(424:478);

% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth2EDC, acc_edc );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure(20);
h = plot( fitresult, xData, yData );
legend( h, 'acc_edc vs. depth2EDC', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth2EDC
ylabel acc_edc
grid on; close figure 20;
ACC_EDC=fitresult(depth_edc); flux_edc=Be10_edc.*ACC_EDC*10^3;
%%
acc_edml=data2.data.EDML(:,6); acc_edml=acc_edml.*(d/an); acc_edml=acc_edml(856:885);
yr_edml=data2.data.EDML(:,2)-1950; yr_edml=yr_edml(856:885);
depth2EDML=data2.data.EDML(:,1); depth2EDML=depth2EDML(856:885);
% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth2EDML, acc_edml );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure(21);
h = plot( fitresult, xData, yData );
legend( h, 'acc_edml vs. depth2EDML', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth2EDML
ylabel acc_edml
grid on; close figure 21
ACC_EDML=fitresult(depth_edml); flux_edml=Be10_edml.*ACC_EDML*10^3;

%%
acc_vostok=data2.data.Vostok(:,6); acc_vostok=acc_vostok.*(d/an); acc_vostok=acc_vostok(291:311);
yr_vostok=data2.data.Vostok(:,2)-1950; yr_vostok=yr_vostok(291:311);
depth2Vostok=data2.data.Vostok(:,1); depth2Vostok=depth2Vostok(291:311);
% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( depth2Vostok, acc_vostok );

% Set up fittype and options.
ft = 'linearinterp';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure(23);
h = plot( fitresult, xData, yData );
legend( h, 'acc_vostok vs. depth2Vostok', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel depth2Vostok
ylabel acc_vostok
grid on; close figure 23
ACC_Vostok=fitresult(depth_vostok); flux_vostok=Be10_vostok.*ACC_Vostok*10^3;
%% Figures of 10Be fuxes
figure(7);
plot(fit1,flux_ngrip,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be flux [atoms/m^2*s]');
hold on; %window_size = 12;
flux_ngrip(isnan(flux_ngrip))=0;
modif = tsmovavg(flux_ngrip,'m',window_size,1);  
plot(fit1,modif,'k', 'Linewidth',1); hold on; 
grid on; box on;
   %title('NGRIP'); %
   title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off

figure(8);
plot(fit2,flux_vostok,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be flux [atoms/m^2 s]');
hold on; %window_size = 12;
flux_vostok(isnan(flux_vostok))=0;
modif22 = tsmovavg(flux_vostok,'m',window_size,1);  
plot(fit2,modif22,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('Vostok,GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off

figure(9);
plot(fit3,flux_edc,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be flux [atoms/m^2 s]');
hold on; %window_size = 12;
flux_edc(isnan(flux_edc))=0;
modif33 = tsmovavg(flux_edc,'m',window_size,1);  
plot(fit3,modif33,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('EDC, GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off

figure(10);
plot(fit4,flux_edml,'b'); hold on; xlabel('year b2k'); 
ylabel('^{10}Be flux [atoms/m^2 s]');
hold on; %window_size = 12;
flux_edml(isnan(flux_edml))=0;
modif44 = tsmovavg(flux_edml,'m',window_size,1);  
plot(fit4,modif44,'k', 'Linewidth',1); hold on; 
grid on; box on;
   title('EDML, GICC05 timescale'); %title('NGRIP,GICC05 timescale')
   legend('data','moving average')
      h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off


figure(11); hold on;
% plot(fit1,flux_ngrip,'b', 'Linewidth',1); hold on;
% plot(fit2,flux_vostok,'r', 'Linewidth',1); hold on;
% plot(fit3,flux_edc,'g', 'Linewidth',1); hold on;
% plot(fit4,flux_edml,'y', 'Linewidth',1); hold on;
plot(fit1,flux_ngrip,'b'); hold on;
plot(fit2,flux_vostok,'r'); hold on;
plot(fit3,flux_edc,'g'); hold on;
plot(fit4,flux_edml,'k'); hold on;
grid on; box on; hold on;  xlabel('year b2k'); ylabel('^{10}Be flux [atoms/m^2 s]'); hold on;
legend('NGRIP (Greenland)', 'Vostok (Antarctica)', 'EDC (Antarctica)', 'EDML (Antarctica)'); hold on
title(' Time series of ^{10} Be flux'); hold on
axis([3.75*10^4 4.55*10^4 0 300]); hold on;       h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',13);
hold off;