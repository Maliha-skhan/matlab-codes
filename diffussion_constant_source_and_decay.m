clear all
close all
clc
n = 5000;
N = 128;
L = 2*pi;
h = L/(N-1);
u(1:N) = 0;
u(floor(N/2)) = 1/h;
u=u' ;
d(1:N)=0 ; %dirac delta function
d(floor(N/2))=1/h ;
d=d' ;
%d(floor(N/2))=1 %dirac delta function
D = 0.1;
t = 0;
tau=1;
S=1;
%dt = (h^2/(2*D))*0.01;
dt=(1/((2*D/(h^2))+1/tau))*0.5
r = D*dt/(h*h); % Stability parameter (r=<0.5)
x0=h*(N/2);
for j = 1:n
   time(1+j)=j*dt;
   
    for i = 1:N
        x(i)=i*h;
        e_ana(i)=exp((-(x(i)-x0)^2)/(sqrt(tau*D)));
        e_ana(i)=abs(e_ana(i));
        if i == 1
           u(i,j+1)=u(i,j)+r*(u(i+1,j)-2*u(i,j))+dt*S*d(i)-(dt/tau)*u(i,j);
           
        elseif i == N
            u(i,j+1)=u(i,j)+r*(u(i-1,j)-2*u(i,j))+dt*S*d(i)-(dt/tau)*u(i,j);

        else
           u(i,j+1) =u(i,j)+r*(u(i-1,j)+u(i+1,j)-2*u(i,j))+dt*S*d(i)-(dt/tau)*u(i,j); % with constant source and decay
           
        end
    end
end
for i=1:N
    e(i)=u(i,n/4)/u(N/2,n/4);
    %e(i)=abs(e(i));
end
%x1=[-N/2 N/2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
plot(x,e_ana,'b')
plot(x,e,'r')
xlabel('x position')
ylabel('u/u_0 (ratio of the temperature at the position of the peak)')
title('Comparison of the solution with the analytical solution at steady state')
legend('Analytical solution','Numerical solution')
box on 
grid on
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%