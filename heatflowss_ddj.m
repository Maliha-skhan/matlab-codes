% Program to calculate steady state temperature profiles 
%
H=3000; % [m]
a=0.06; % [m/yr]
Qgeo=35/1000; % [W/m2]
Tsur=-30;  %{deg C]
year2second=365.25*24*60*60;
N=51;
%
% Constants
C=2097;  %[Jkg-1K-1]
L=333.5*1000;  %[Jkg-1]
K=2.10;   %[Wm-1K-1]
rho_ice = 918; %[kg/m3]
kappa=K/(rho_ice*C);  % [m2s-1]
%
del_z=H/(N-1);
z=[0:del_z:H]';
A=zeros(N,N);
B=zeros(N,1);
T=zeros(N,1);
w=-([1:N]'-1)*a/(N-1); %[m/yr]
w=w/year2second;    %[m/s]

%
% Define A and B - assuming Tbase > Tmelt
for i=2:N-1
    A(i,i)=-2;
    A(i,i-1)=1+del_z*w(i)/(2*kappa);
    A(i,i+1)=1-del_z*w(i)/(2*kappa);
end
% Boundary condition, base
A(1,1)=1;
A(1,2)=-1;
B(1)=del_z*Qgeo/K;
% Boundary condition, surface
A(N,N)=1;
B(N)=Tsur;
%
% Invert AT=B
T=inv(A)*B;  %[or T=(A^(-1))*B;]
figure(1)
plot(T,z,'b')
hold on
grid on
%
% change to melt conditions if the basal temperature is above 0
if T(1)>0
    A(1,1)=1;
    A(1,2)=0;
    B(1)=0;
    % Invert AT=B
    T1=inv(A)*B;  %[or T=(A^(-1))*B;]
    Qice=-K*(T1(2)-T1(1))/del_z;
    Qmelt=Qgeo-Qice;
    wmelt=-Qmelt/(rho_ice*L);
    wmelt_myr=-wmelt*year2second;
    plot(T1,z,'r')
    %adjust w
    w1=-([1:N]'-1)*(a+wmelt_myr)/(N-1)+wmelt_myr; %[m/yr]
    w1=w1/year2second;    %[m/s]
    for i=2:N-1
    A(i,i)=-2;
    A(i,i-1)=1+del_z*w1(i)/(2*kappa);
    A(i,i+1)=1-del_z*w1(i)/(2*kappa);
    end
    % Invert AT=B
    T2=inv(A)*B;  %[or T=(A^(-1))*B;]
    plot(T2,z,'g')
end
    
    
    

