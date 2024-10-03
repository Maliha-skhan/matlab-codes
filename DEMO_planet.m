clear all
close all
clc

E_M = 1 ; % earth mass
Sun_M=1000 ; % sun mass
G=1*10^-3 ;
% Earth
x1=[5.5 0]
v1=[0 0.4]

%Mars
m2=0.9
v2=[0 0.23]
x2=[9.2 0.1 ]

% Mercury
m3=0.553
x3=[3.15 0]
v3=[0 0.45]

%Comet
m4=8
x4=[19 0]
v4=[0 0.167895201]

figure(1)
for i=1:700 ;
% sun
plot(0,0,'r*')
    %Earth
    hold on
Dt=0.5 ;
F1=((-G*Sun_M*E_M)/(((norm(x1))^2)))*(x1/norm(x1)) ;
a1=F1/E_M ;
v1=v1+a1*Dt ;
x1=x1+v1*Dt ;

x_p_1(i)=x1(1) ;
y_p_1(i)=x1(2) ;

plot(x_p_1(1:i),y_p_1(1:i),'g-')
plot(x1(:,1),x1(:,2),'g*')

axis([-15 20 -15 15])
% Mars
F2=((-G*Sun_M*m2)/(((norm(x2))^2)))*(x2/norm(x2)) ;
a2=F2/m2 ;
v2=v2+a2*Dt ;
x2=x2+v2*Dt ;

x_p_2(i)=x2(1) ;
y_p_2(i)=x2(2) ;


plot(x_p_2(1:i),y_p_2(1:i),'k-')
plot(x2(:,1),x2(:,2),'k*')

% Mercury
F3=((-G*Sun_M*m3)/(((norm(x3))^2)))*(x3/norm(x3)) ;
a3=F3/m3 ;
v3=v3+a3*Dt ;
x3=x3+v3*Dt ;

x_p_3(i)=x3(1) ;
y_p_3(i)=x3(2) ;

plot(x_p_3(1:i),y_p_3(1:i),'b-')
plot(x3(:,1),x3(:,2),'b*')
title('Solar System: Planetary motion')
xlabel('x-position [random unit]')
ylabel('y-position [random unit]')

% Comet
F4=((-G*Sun_M*m4)/(((norm(x4))^2)))*(x4/norm(x4)) ;
a4=F4/m4 ;
v4=v4+a4*Dt ;
x4=x4+v4*Dt ;

x_p_4(i)=x4(1) ;
y_p_4(i)=x4(2) ;

if i<50
plot(x_p_4(1:i),y_p_4(1:i),'m-')
plot(x4(:,1),x4(:,2),'m*')
end
if i>49
plot(x_p_4(i-49:i),y_p_4(i-49:i),'m-')
plot(x4(:,1),x4(:,2),'m*')
end
hold off
pause(0.01)
title('Solar System: Planetary motion')
xlabel('x-position [random unit]')
ylabel('y-position [random unit]')
%legend('Sun','Earth', 'Mars', 'Mercury', 'location', 'northwest')
end
figure(2)
plot(0,0,'r*')
hold on
plot(x_p_1,y_p_1,'g.')
plot(x_p_2,y_p_2,'k.')
plot(x_p_3,y_p_3,'b.')
plot(x_p_4,y_p_4,'m.')

axis([-15 20 -15 15])
title('Solar System: Planetary motion')
xlabel('x-position [random unit]')
ylabel('y-position [random unit]')
legend('Sun','Earth', 'Mars', 'Mercury', 'Comet', 'location', 'northwest')
