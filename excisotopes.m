clear all; close all; clc;
d18O =[-25.5 -40.5 -28.0] ;
temp=10^(-1)*[-100 -400 -150] ;

dTperdox= (-105+140)/(10*(-25.6+27.5))

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( temp, d18O );

% Set up fittype and options.
ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.0975404049994095 0.278498218867048];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'd18O vs. temp', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel temp
ylabel d18O
grid on

fitresult

%% d) I notesættets figur 5 ses data fra 0-123.000 år b2k. Hvornår var det koldest/varmest?
% koldest ca. 28 000 b2k
% varmest 123 000 b2k




%%%Exercise on changing stable isotope values in the ocean in glacial periods
% Question 1
m_ocean=  1.35*10^9*(1000)^3 * 1000 % kg
m_ice = (2.95+29+0.18+0.05)*10^6*(1000)^3 * 918 % kg
frac=m_ice/m_ocean

% Question 2
d18O_standardocean=0
d18O_ice=-30 % per mille

% Question 3
%weights
w_oceanglacial=3.6*10^8*(1000)^2*(3700-130)* 1000
w_iceglacial=(m_ocean+m_ice)-w_oceanglacial
%change 
ch_ice=m_ice-w_iceglacial
ch_ocean=m_ocean-w_oceanglacial

% Question 4
%Assume the mean ?18O value of the ice sheets are unchanged during the glacial period compared to present
%Use the ratio of 18O to 16O to make a balance between the present and the glacial content of 18O
ratio_present=0.9700 % present value for ([18O]/[16O])

% Question 5
% What would the mean ?18O of the ocean be during the glacial? 
ratio_glacial=ratio_present*m_ocean*(w_oceanglacial)^(-1)
(ratio_glacial-1)*1000

% 6. Compare the result with the curve from the article of Bintanja in Nature 2005.  
 