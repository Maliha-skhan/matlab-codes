clear all;close all; clc;
%%%Exercise: The Vertical Fault { Part 2
%%%Linear inversion through Tikhonov regularization
rho_homo = 2600;
G=6.67*10^(-11); % is the Gravity constant
gravc=G
% load data
gravdata=load('gravdata.txt');
dj=gravdata(:,2); % this corresponds to our d matrix
xj=1000*gravdata(:,1); %  converted into km

% In this assignment we represent the subsurface to the right of the fault by only 5 horizontal layers.
% This time, not only the mass density contrasts, but also the thicknesses of the layers are unknown.

%d_j=G*rho_dif*log((zbase^2+x^2)/(ztop^2+x^2))
z=0:10^3:10^3*10; eps=0.05;
%x=0:1:18;
for i=1:length(z);
    for j=1:length(xj);
%   g(i,j) = log10((z(i)^2 + x(j))/( ((z(i)+1)^2 + x(j) ))); % this is our G matrix
% g(j,i) = G * log10((z(i)^2 + xj(j))/( ((z(i)-1000)^2 + xj(j) ))) ...
%     + dgdm(j,i)*dm ; % this is our G matrix
gold(j,i) = G * log10((z(i)^2 + xj(j))/( ((z(i)-1000)^2 + xj(j) ))); % this is our G matrix

    end
end
m_old = ( transpose(gold)*gold + eps^2*eye(length(z)) )^(-1)*transpose(gold)*dj
m0=[45e-14 41e-14 24e-14 11e-14 6e-14 0 3 5 8 13]' % my guesses


G = zeros(18,5);
m_rho=m0(1:5);
m_z=zeros(6,1);
m_z(2:6)=m0(6:10);
for j=2:6
    for i=1:18;
        G(i,j-1)= gravc*log10( (m_z(j)^2 + xj(i)^2)/(m_z(j-1)^2 + xj(i)^2 ) );
    end
end
d0=G*m_rho;
 g_m(m0,xj)
% the function if working correctectly!!

dm0=[0.00001*45e-14 0.00001*41e-14 0.00001*24e-14 0.00001*11e-14 0.00001*6e-14 0.00001 0.00001*3 0.00001*5 0.00001*8 0.00001*13]' % my guesses
for i=1:10
    for j=1:10
        dm(i,i)=dm0(i)
dgdm(i)= (g_m(m0+dm,xj)-g_m(m0,xj))/(dm(i,i));
    end; end

% path
% which g_m.m -all

%d0=g_m(m_old,xj)


% for k=1:10;
%     dm(k)= 0.001*ones(1,10) %m_mold
%     dgdm(k)=(gold*(m_old - dm(k))-gold*(dm_old))/dm(k)
% end

% Compute the resolution matrix for the vertical fault problem, and explain the result.
% resolution matrix R = [GTG + 2I]?1GTG,
R_old =( transpose(gold)*gold + eps^2*eye(length(z)) )^(-1)*transpose(gold)*gold;
% resolution matrix tells us 
%The resolution matrix is, so to speak, the "window" through which we see the true model mtrue. 
%The "picture" m_e is usually blurred, and the columns of R are the "blurring functions" describing 
%how the details (isolated `spikes') are smeared out in the Damped Least Squares solution




% this is a under-determined problem so it shouldn't be unique

% 
% figure(1); plot(z,m_e,'mx','linewidth',2); xlabel('depth [km]'); ylabel('\Delta\rho'); hold on;
% title('the vertical density variation');box on; grid on; %hold on
% 
% %plot(z,R(:,10),'g-'); legend('density','resolution column 10')
% 
% figure(2); plot(z,R(:,10),'g-','linewidth',2);hold on; title('resolution column 10');  
% xlabel('depth [km]'); ylabel('resolution');box on; grid on;


% function mest = gmatfun(x,z,d)
%  G=10;
%  for i=1:length(z);
%     for j=1:length(x);
% %   g(i,j) = log10((z(i)^2 + x(j))/( ((z(i)+1)^2 + x(j) ))); % this is our G matrix
% % g(j,i) = G * log10((z(i)^2 + xj(j))/( ((z(i)-1000)^2 + xj(j) ))) ...
% %     + dgdm(j,i)*dm ; % this is our G matrix
% gold(j,i) = G * log10((z(i)^2 + x(j))/( ((z(i)-1000)^2 + x(j) ))); % this is our G matrix
% %dgdm(j,i) = 
%     end
%  end
% mest = ( transpose(gold)*gold + eps^2*eye(length(z)) )^(-1)*transpose(gold)*d;
%  end