function [H,Lh] = cons_hypergraph(X,k)
m = size(X,1);
H = zeros(m,m);    
for j = 1:m              
    z = knnsearch(X,X(j,:),'K',k);  
    for i = 1:k
         H(z(i),j) = 1;
    end
end
A = eye(m);    
Dv = zeros(m,m);   
tmp1 = sum(H'*A,2);
for i = 1:m
    Dv(i,i) = tmp1(i);
end
De = k*eye(m);  
theta1 = (Dv^-0.5)*H*A*(De^-1)*H'*(Dv^-0.5);
I = eye(m);
Lh = I - theta1;    % LapLacian matrix;