clear all; close all; clc
%-----------------------------------------------------------------------%
%-----------------The Quantum Harmonic Oscillator-----------------------%
%-----------------------------------------------------------------------%


N=400;  % no of grid points
x=linspace(-20,20,N); % x-coordinate % -L : L 

psi0=exp(-x.*x); % ground state wave function, Griffiths, p. 58, Eq. 2.59
%%% While making all constants equals 1

%---forming the raising operator ------------------------------------%

neg_1=-1.*ones(N,1); % Only containing -1 
pos_1=1.*ones(N,1);  % Only containing +1 
X=spdiags(x,N,N); % x operator
P=0.5.*(spdiags([neg_1 pos_1 ],[-1 1] ,N,N)); %momentum operator

aplus=(sqrt(2))^-1.*(X-P); % Griffiths, p. 54, Eq. 2.47


psi0=psi0-mean(psi0);  
psi0=psi0./max(psi0);
psi(:,1)=psi0;
%--- Using the raising operator to find excited states------------------%
for i=2:40
    psi(:,i)=(aplus).^(i) * psi(:,i-1);
    psi(:,i)=psi(:,i) - mean(psi(:,i)); % normalizing
    psi(:,i)=(psi(:,i)./max(psi(:,i)));
end

V=0.5*x.*x; % parabolic oscillator potential
V=V./max(V);
V=100.*V;

%----------------------------------------------------------------------%

figure(1)
%------rescaling the wavefunction for visualization purpose----------%
plot(x,psi(:,1),'b','linewidth',2)
   text(0, 1.5,'E_0','fontSize',14)
axis([-5 5 0 15])
hold on
plot(x,psi(:,2)+3,'b','linewidth',2)
   line([-10 10], [3 3],'color','k','Linestyle','-.')
   text(-0.5, 3.5,'E_1','fontSize',14)
plot(x,psi(:,3)+7,'b','linewidth',2)
   line([-10 10], [7 7],'color','k','Linestyle','-.')
   text(0.5, 7.5,'E_2','fontSize',14)
plot(x,psi(:,4)+12,'b','linewidth',2)
   line([-10 10], [12 12],'color','k','Linestyle','-.')
   text(0, 12.5,'E_3','fontSize',14)
plot(x,V,'r','linewidth',2)

   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('\psi','fontSize',14);
   title('Stationary States of Harmonic Oscillator','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   set(gca,'ytick',[]);

hold off
figure(2)
%------rescaling the wavefunction for visualization purpose----------%
plot(x,psi(:,10),'b','linewidth',2)
   text(0, 1.5,'E_{11}','fontSize',14)
%axis([-5 5 0 15])
hold on
plot(x,psi(:,11),'b','linewidth',2)
   line([-10 10], [3 3],'color','k','Linestyle','-.')
   text(-0.5, 3.5,'E_{12}','fontSize',14)
plot(x,psi(:,12),'b','linewidth',2)
   line([-10 10], [7 7],'color','k','Linestyle','-.')
   text(0.5, 7.5,'E_{13}','fontSize',14)
plot(x,psi(:,13),'b','linewidth',2)
   line([-10 10], [12 12],'color','k','Linestyle','-.')
   text(0, 12.5,'E_{14}','fontSize',14)
plot(x,V,'r','linewidth',2)

   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('\psi','fontSize',14);
   title('Stationary States of Harmonic Oscillator','fontsize',14)
   fh = figure(1);
   set(fh, 'color', 'white'); 
   set(gca,'ytick',[]);

hold off
%------------------------end--------------------------------------------%
%------------------------end--------------------------------------------%