clear all
close all
clc

E_M = 1
G=1*10^-3
E_R=1
Sun_M=1000
E_F=G*Sun_M*E_M/((E_R)^2)
a_e=E_F/E_M
w1=-sqrt(a_e/E_R)
t=0
x=0

figure(1)
hold


%w1=1

for i=1:200 ;
t(i)=(i-1)*2*pi/100;
%planet 1
x1(i)=E_R*cos(w1*t(i))
y1(i)=E_R*sin(w1*t(i))
plot(x1,y1,'ok')
%planet 2
x2(i)=0.4*E_R*cos(2*w1*t(i))
y2(i)=0.4*E_R*sin(2*w1*t(i))
plot(x2,y2,'or')

%planet 3
x3(i)=1.8*E_R*cos(0.6*w1*t(i))
y3(i)=1.8*E_R*sin(0.6*w1*t(i))
plot(x3,y3,'or')
axis([-3 3 -3 3])
%comet(x,y)
plot(0,0,'o')
pause(0.05)
%axis([0 10^16 0 10^16])
end
