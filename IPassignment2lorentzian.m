clear all; clc; close all;
load mars_soil.txt
%[z d] = [mars_soil(:,1) mars_soil(:,2)];
z=mars_soil(:,1); d =mars_soil(:,2); d =max(d)-d;
plot(z,d,'r');hold on; title('Data: mars soil'); grid on; 
box on;xlabel('velocity [mm/s]'); ylabel('counts')

% overedetermined problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The Lorentzian case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the g(m) matrix
%dA,df, dc
coef=1/sqrt(2*pi);
%A1=100; f1=0.0147; c1=2.650; % coefficients going inside matrix m
% constructing matrix m for p=1
pp=17;
m=zeros(3*pp,1);
m(1,1)=.5*500; 
m(2,1)=-11; 
m(3,1)=0.9; 
m(4,1)=.5*1001;
m(5,1)=-8.8;
m(6,1)=0.80;
m(7,1)=.5*1004.4;
m(8,1)=-7.6;
m(9,1)=0.2;
m(10,1)=.5*905.5;
m(11,1)=-6.521;
m(12,1)=0.07;
m(13,1)=.5*1222;
m(14,1)=-5.5;
m(15,1)=0.05;
m(16,1)=.5*50; 
m(17,1)=-4.3; 
m(18,1)=0.1; 
m(19,1)=.5*821;
m(20,1)=-3.7;
m(21,1)=0.030;

m(22,1)=.5*720;
m(23,1)=-3.1;
m(24,1)=0.580;
m(25,1)=.5*524.4;
m(26,1)=-1.6;
m(27,1)=0.12;
m(28,1)=.5*905.5;
m(29,1)=1.721;
m(30,1)=0.07;
m(31,1)=.5*822;
m(32,1)=3.2;
m(33,1)=0.05;
m(34,1)=.5*770; 
m(35,1)=4.3; 
m(36,1)=0.9; 
m(37,1)=.5*720;
m(38,1)=5.4;
m(39,1)=0.80;
m(40,1)=.5*624.4;
m(41,1)=6.4;
m(42,1)=0.2;
m(43,1)=.5*965.5;
m(44,1)=7.521;
m(45,1)=0.07;
m(46,1)=.5*1132;
m(47,1)=8.7;
m(48,1)=0.05;
m(49,1)=.5*132;
m(50,1)=10.5;
m(51,1)=0.5;




%[A1 f1 c1];
%g(m) for lorentzian
%glos(:,1) = ( m(1)*m(3)^2 )/(( z- m(2) ).^2+m(3)^2);
% for i=1:length(z)
% Gmatrix(i,1) = m(3,1)^2/( (z(i)-m(2,1))^2 +m(3,1)^2); %ok
% Gmatrix(i,2) = (2*m(1,1)*(z(i)-m(2,1)) *m(3,1)^2 )/( (z(i)-m(2,1))^2 +m(3,1)^2)^2;
% Gmatrix(i,3) = (2*m(1,1)* m(3,1) *(z(i)-m(2,1))^2 )/( (z(i)-m(2,1))^2 +m(3,1)^2)^2;
% glos(i,1) = ( m(1,1)*m(3,1)^2 )/(( z(i)- m(2,1) )^2+m(3,1)^2);
% end

%alpha=0.0005; %  between 0 and 1, but can be made time-dependent

alpha=0.001; 
% delta_m(1,1)=0;
% delta_m(2,1)=0;
% delta_m(3,1)=0;
delta_m=zeros(pp*3,1);

for n=1:1000
for i=1:length(z);
   for p=1:pp;
%     glo(i) =  ( m(1,1)*m(1,3)^2 )/(( z(i)- m(1,2) )^2+m(1,3)^2);
%     Gmatrix(i,1) = m(1,3)^2/( (z(i)-m(1,2))^2 +m(1,3)^2);
% Gmatrix(i,2) = (2*m(1,1)* m(1,3)^2 *(z(i)-m(1,2)) )/( (z(1)-m(1,2))^2 +m(1,3)^2)^2;
% Gmatrix(i,3) = (2*m(1,1)* m(1,3) *(-z(i)+m(1,2))^2 )/(m(1,3)^2 + m(1,2)^2-z(i)*m(1,2)+ z(1)^2);

% m(3*p-2)= m(3*(p-1)-2)+alpha*delta_m(3*(p-1)-2);
% m(3*p-1)= m(3*(p-1)-1)+alpha*delta_m(3*(p-1)-1);
% m(3*p)= m(3*(p-1))+alpha*delta_m(3*(p-1)); %[A1 f1 c1];

% m(3*p-2,1)= m(3*(p-1)-2,1)+alpha*delta_m(3*p-2);
% m(3*p-1,1)= m(3*(p-1)-1,1)+alpha*delta_m(3*p-1);
% m(3*p,1)= m(3*(p-1),1)+alpha*delta_m(3*p); %[A1 f1 c1];

%glos(i,p) =  ( m(3*p-2,1)*m(3*p,1)^2 )/(( z(i)- m(3*p-1,1) )^2+m(3*p,1)^2);
glos(i,p) =  ( m(3*p-2,1)*m(3*p,1)^2 )/(( z(i)- m(3*p-1,1) )^2+m(3*p,1)^2);


% glos(i,1) = ( m(1,1)*m(3,1)^2 )/(( z(i)- m(2,1) )^2+m(3,1)^2);

Gmatrix(i,3*p-2) = m(3*p,1)^2/( (z(i)-m(3*p-1,1))^2 +m(3*p,1)^2);
Gmatrix(i,3*p-1) = (2*m(3*p-2,1)* m(3*p,1)^2 *(z(i)-m(3*p-1,1)) )/( (z(i)-m(3*p-1,1))^2 +m(3*p,1)^2)^2;
Gmatrix(i,3*p) = (2*m(3*p-2,1)* m(3*p,1) *(z(i)-m(3*p-1,1))^2 )/( (z(i)-m(3*p-1,1))^2 +m(3*p,1)^2)^2;



   end
   if i==length(z);
       glo=sum(glos,2);
       glo=glo';
       delta_m = (transpose(Gmatrix)*Gmatrix)^(-1)*transpose(Gmatrix)*(d-glo);
   end
end
for ii=1:pp*3
m(ii) = m(ii)+alpha*delta_m(ii); %finish later
end
end
%    if i==length(z);
% glo=sum(glos,2);
% glo=glo';
% delta_m = (transpose(Gmatrix)*Gmatrix)^(-1)*transpose(Gmatrix)*(d-glo);

% for n=1:10
%     
% end

%    end

% glo=glo';
% delta_m = (transpose(Gmatrix)*Gmatrix)^(-1)*transpose(Gmatrix)*(d-glo);

% for p=2:20
%     m(p,:)= m(p-1,:)+alpha*delta_m;
% end
 figure();hold on; plot(glo,'g','Linewidth',1);hold on;plot(d,'k');