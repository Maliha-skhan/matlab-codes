clear all; close all; clc;
wavedata=load('wavedata.txt');
pointspread=load('pointspread.txt');
P=pointspread; D=wavedata;
for i=1:50
for j=1:50
    D_matrix(i,j)=D(i,j);
    d(i*j,1)=D_matrix(i,j);
end
end
id=eye(100);
%k=rand(400,1);
k=ones(400,1); 
M=reshape(k,[20,20]);
Deltamodel=zeros(20,20); Deltamodel(1,4)=1; Deltamodel(5,6)=1; Deltamodel(12,18)=1; 
Deltamodel(17,20)=1; Deltamodel(15,5)=1;
% Dk = conv2(M,P,'same');
% G=Dk
Dk = conv2(Deltamodel,P,'same');
G=Dk;
gk=reshape(Dk, [(size(Dk,1)*size(Dk,1)),1]);

% Find a solution m_e to the problem using Tikhonov Regularization
 eps=0.2;
%eps=sqrt(300);

m_e=zeros(size(d,1),1);

for n=1:size(d,1)/20
m_e(1+20*(n-1):20+20*(n-1)) = ( transpose(G)*G + eps^2*eye(20) )^(-1)*transpose(G)*d(1+20*(n-1):20+20*(n-1));
end

mm=reshape(m_e, [50,50]);
figure();imagesc(mm)

% Compute the resolution matrix for the vertical fault problem, and explain the result.
% resolution matrix R = [GTG + 2I]?1GTG,
%R=( transpose(P)*P + eps^2*eye(101) )^(-1)*transpose(P)*P 

%plot(dk)
%eps=0.0005
%m = ( transpose(P)*P + eps^2*id )^(-1)*transpose(P)*d

%M=((D^(-1))*P)^(-1)
% for k=1: size(D,1);
%     
%     M()
% end
% figure(); mesh(P); title('Point spread')
% figure(); mesh(Dk); title('Dk')
% figure(); mesh(D); title('untouched wave data')