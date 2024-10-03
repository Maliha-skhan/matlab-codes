clear all;close all; clc;
%%%Exercise: The Vertical Fault { Part 1
%%%Slightly nonlinear inversion through Tikhonov regularization
rho_homo = 2600;
G=6.67*10^(-11); % is the Gravity constant

% load data
gravdata=load('gravdata.txt');
dj=gravdata(:,2); % this corresponds to our d matrix
xj=1000*gravdata(:,1); % 

%d_j=G*rho_dif*log((zbase^2+x^2)/(ztop^2+x^2))
z=0:10^3:10^3*99;
x=0:1:18;
for i=1:100;
    for j=1:18;
%   g(i,j) = log10((z(i)^2 + x(j))/( ((z(i)+1)^2 + x(j) ))); % this is our G matrix
g(j,i) = G * log10((z(i)^2 + xj(j))/( ((z(i)-1000)^2 + xj(j) ))); % this is our G matrix
    end
end

%g_sum=G*sum(g,1);


sigma=10^(-9); %Data is measured with an uncertainty sigma


% Find a solution m_e to the problem using Tikhonov Regularization
eps=0.2


m_e = ( transpose(g)*g + eps^2*eye(100) )^(-1)*transpose(g)*dj

% Compute the resolution matrix for the vertical fault problem, and explain the result.
% resolution matrix R = [GTG + 2I]?1GTG,
R=( transpose(g)*g + eps^2*eye(100) )^(-1)*transpose(g)*g 
% resolution matrix tells us 
%The resolution matrix is, so to speak, the "window" through which we see the true model mtrue. 
%The "picture" m_e is usually blurred, and the columns of R are the "blurring functions" describing 
%how the details (isolated `spikes') are smeared out in the Damped Least Squares solution




% this is a under-determined problem so it shouldn't be unique


figure(1); plot(z,m_e,'mx','linewidth',2); xlabel('depth [km]'); ylabel('\Delta\rho'); hold on;
title('the vertical density variation');box on; grid on; %hold on

%plot(z,R(:,10),'g-'); legend('density','resolution column 10')

figure(2); plot(z,R(:,10),'g-','linewidth',2);hold on; title('resolution column 10');  
xlabel('depth [km]'); ylabel('resolution');box on; grid on;