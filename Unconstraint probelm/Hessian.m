function H= Hessian(A,x)
[m,n]=size(A);
H=zeros(n,n);
for i=1:m
    H=H+A(i,:)'*A(i,:)/((1-A(i,:)*x)^2);
end
for j=1:n
    H(j,j)=H(j,j)+(1/((1+x(j,1))^2)+1/((1-x(j,1))^2));
end