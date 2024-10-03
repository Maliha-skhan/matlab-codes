close all
clear all
clc

%Moon
G=6.7*10^-11
E_M=6*10^24
Dt=0.1
m4=1
x4=[5.15 0]
v4=[0 0.40002679]
R_em=[7 0]
figure(3)
for i=1:500
% Moon
hold on
F4=-G*m4*E_M/((norm(R_em)^2))*(R_em/norm(R_em)) ;
a4=F4/m4 ;
v4=v4+a4*Dt ;
%x4=x4+v4*Dt ;
%R_em=x4-x1
R_em=R_em+v4*Dt
%F4=((-G*Sun_M*m4)/(((norm(x4))^2)))*(x4/norm(x4))+((-G*m4*E_M)/(((norm(R_em))^2)))*(-R_em/norm(R_em)) ;

x_p_4(i)=R_em(1) ;
y_p_4(i)=R_em(2) ;
xR(i)=R_em(1) ;
yR(i)=R_em(2) ;
pause(0.01)
%plot(x_p_4(1:i),y_p_4(1:i),'m-')
%plot(x4(:,1),x4(:,2),'m*')
plot(xR(1:i),yR(1:i),'m-')
plot(R_em(:,1),R_em(:,2),'m.')
axis([-15 15 -15 15])
end