clear all; close all; clc;
dxdt1=[2/1 2/1.5 2/2];
error1=[1 4 5];

dxdt2=[2/2 5/2 10/2];
error2=[5 16 33];

dxdt31=[2/1.6667 10/1.6667 20/1.6667 23/1.6667];
error31=[1 19 40  287];

dxdt3=[2/1.6667 10/1.6667 20/1.6667];
error3=[1 19 40];

dxdt=[dxdt1 dxdt2 dxdt3];
error=[error1 error2 error3];

dxdt11=[dxdt1 dxdt2 dxdt31];
error11=[error1 error2 error31];

figure(1)
plot(dxdt,error,'bo'); grid on; box on
xlabel('\Delta x/ \Delta t [km/min]')
ylabel('Percentage error [%]')
title('Stability of the numerical model')

figure(3)
 %Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( dxdt11, error11 );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.756081594811466 0.383063807777644 0.2];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'error11 vs. dxdt11', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel dxdt11
ylabel error11
grid on
%%
figure(22)
plot(dxdt11,error11,'bo'); hold on; grid on; box on
plot(fitresult); legend('Data point','Fitted line'); hold on
text(5,250,'y_{fit}=9.90*10^{-6} * exp(1.24x) + 10.49','HorizontalAlignment','left', 'fontsize', 12);
xlabel('\Delta x/ \Delta t [km/min]')
ylabel('Percentage error [%]')
title('Stability of the numerical model');
set(gca,'fontsize', 12);set(gca,'DefaultTextFontSize',12); 
