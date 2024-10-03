clear all
close all
clc
n = 80000;
N = 128;
L = 2*pi;
h = L/(N-1);
u(1:N) = 0;
u(floor(N/2)) = 1/h;
d(1:N)=0 %dirac delta funfuktion
%d(floor(N/2))=1/h
d(floor(N/2))=1 %dirac delta funfuktion
D = 0.1;
t = 0;
tau=1;
S=1;
%dt = (h^2/(2*D))*0.001;
dt=0.5*(2/((4*D/(h^2))+(1/(tau))));
u_ana(1:N) = 0;
u_ana(N/2) = 1/h;
r = D*dt/(h*h); % Stability parameter (r=<0.5)
u=u';
u_ana = u_ana';
d=d'
x0=h*N/2
for j = 1:n
   time(1+j)=j*dt;
        sigma(j+1) = sqrt(2*D*j*dt);
    for i = 1:N
        u_ana(i,j+1) = 1/((sigma(j+1))*sqrt(2*pi))*exp(-((h*(i-N/2))^2)/(2*sigma(j+1)^2));
        x(i)=i*h;
      
        e_ana(i)=exp(-abs((i*h)-x0)/(sqrt(tau*D)));
        if i == 1
           u(i,j+1)=u(i,j)+r*(u(i+1,j)-2*u(i,j))-(dt/tau)*u(i,j)+dt*S*d(i);
           %e(i,j+1)=u(i,j+1)/u(1,j+1);
           
        elseif i == N
            u(i,j+1)=u(i,j)+r*(u(i-1,j)-2*u(i,j))-(dt/tau)*u(i,j)+dt*S*d(i);
          %  e(i,j+1)=u(i,j+1)/u(1,j+1);
        else
           u(i,j+1) =u(i,j)+r*(u(i-1,j)+u(i+1,j)-2*u(i,j))-(dt/tau)*u(i,j)+dt*S*d(i); % with constant source and decay
         %e(i,j+1)=u(i,j+1)/u(1,j+1);
        end
    end
end
for i=1:N
    e(i)=u(i,n)/u(N/2,n);
    e(i)=abs(e(i));
end

    
Err = abs(u_ana-u);
floor(n/2)
figure(1) 
hold on
plot( u(:,1),'.-r' )
hold on
plot(u_ana(:,1),'.-b')
hold on
plot(u_ana(:,n),'.-g')
hold on
 plot(u(:,n),'.-k')
hold on
 plot(u_ana(:,floor(n/2)),'.-m')
 hold on
plot(u(:,floor(n/2)),'.-c')
hold on
title('Temperature at different times')%  plot(c_ana(:,150))
xlabel('number of grid point in space')
ylabel('Temperature')
legend('u at time 0 sec','u_{ana} at time 0 sec','u_{ana} n/2=0.31 sec','u n/2=0.31 sec','u_{ana} n=0.61 sec','u n=0.61 sec','Location','best')   %'c-ana N/4','c-ana N/8','c N','c N/2','c N/4','c N/8', 
grid on
box on
%axis([58 70  ])
hold off

figure(2)
hold on
plot(Err(:,1))
plot(Err(:,n))
plot(Err(:,floor(n/2)))
%plot(Err(:,floor(n/4)))
%plot(Err(:,floor(n/8)))
%axis([60 70 0 5])
legend('Err at time 0 sec','Err at time n=0.61 sec','Err at time n/2=0.31')  %,'Err N/2','Err N/4','Err N/8')
grid on
box on
hold off

% *Derive von Neumann stability condition and test it
figure(3)
hold on
mesh(x,time,u')
xlabel('x position')
ylabel('time')
zlabel('Temperature')

figure(4)
hold on
plot(x,e_ana,'xb'); hold on
plot(x,e,'r')
title({'Comparison of the solution with the analytical'; 'solution at steady state (FTCS)'})
legend('Analytical solution','Numerical solution'); xlabel('x position'); ylabel('u/u_0 (ratio of temperature)')
box on 
grid on
hold off