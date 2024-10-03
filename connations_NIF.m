clear all
close all
clc


global node
%Ntot = 500;
r = 60;
Ntot=0.5*r*(r-1);
%p = Ntot/((r*(r-1))/2)
p=0.05
N = Ntot*p

for i = 1:r
node(i).neighbor = single.empty;
end



for  k = 1:N
%     if k <= N/2
    i = randi([1 r]);
    j = randi([1 r]);
    if i ~= j
    if isempty(node(i).neighbor) == 0    
        for g =1:length(node(i).neighbor)
            if node(i).neighbor(g) == j
                break
            elseif g == length(node(i).neighbor)
                 node(i).neighbor(end+1) = j;
                 node(j).neighbor(end+1) = i;              
            end
            
        end
    else
        node(i).neighbor(end+1) = j;
        node(j).neighbor(end+1) = i;
    end
    end

    
end
% der plotes
for i = 1:r
vis(1:r,i) = zeros;
vis(node(i).neighbor,i) = 1;
end


G = graph(vis);
plot(G);

% node(1).neighbor=[2 3 4]
% node(2).neighbor=[3 5]
% node(3).neighbor=[2 6]
% node(4).neighbor=[4 6]
% node(5).neighbor=[2]
% node(6).neighbor=[3 4]
% 
% node(7).neighbor=[8]
% node(8).neighbor=[7 9]
% node(9).neighbor=[8]

%% der farves
alist = [];
for i = 1:length(node)
node(i).color = 0;
end

farve = [];
c = 1;
for i = 1:length(node)
    alist = [];
if node(i).color == 0
    c = c+1;
    farve(end+1) =c; 
    node(i).color = c;
%     alist(end+1) = i;
    cnf(i,c);
end
end


bigest = [];

% der tælles 
for i = 2:max(farve)
    for j= 1:r
    if node(j).color == i;
    bigest(i,end+1) = 1;
    end
    end
    bigest_s(i) = sum(bigest(i,:)); 
end
%%


bigest = [];

for k = 1:r
   p = randi([1 r]);
   nej = 1;
   while isempty(node(p).neighbor) == 1 & nej <= 2*r
       p = randi([1 r]);
       nej = nej + 1;
   end
        for g =1:length(node(p).neighbor) 
            gg = node(p).neighbor(g); 
            for y = 1:length(node(gg).neighbor)
                if node(gg).neighbor(y) == p
                    tag = y; 
                end
            end
            node(gg).neighbor(tag) = []
        end
        node(p).neighbor = [];
       
  % de bliver farvet      
    alist = [];
for i = 1:length(node)
node(i).color = 0;
end

farve = [];
c = 1;
for i = 1:length(node)
    alist = [];
if node(i).color == 0
    c = c+1;
    farve(end+1) =c; 
    node(i).color = c;
%     alist(end+1) = i;
    cnf(i,c);
end
end

% der tælles   
bigest = []; 
bigest_s = [];
for i = 2:max(farve)
    for j= 1:r
    if node(j).color == i;
    bigest(i,end+1) = 1;
    end
    end
    bigest_s(i) = sum(bigest(i,:)); 
end   
clost(k) = max(bigest_s)/r; 


end

RN = (1:r)/r;
figure(2); hold on; plot(RN,clost); ylabel('S_{LC}/N_{tot}'); xlabel('Percentage of removed nodes');
title('The size of the giant component as a function of number of removed nodes/edges'); grid on; box on; hold off 