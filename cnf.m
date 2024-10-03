function  farve = cnf(M,c) 

 global node
alist = []; 
alist(end+1) = M;
while isempty(alist) == 0
    kk = alist(1);
     
    for j = 1:length(node(kk).neighbor)
         jj = node(kk).neighbor(j);
         if node(jj).color == 0
           alist(end+1)=jj;
           node(jj).color = c;     
         end
        
    end
    alist(1) = [];
end

  
 

end