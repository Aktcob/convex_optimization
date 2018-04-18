function grad_f= grad_obj_func(A,x)
[m,n]=size(A);
grad_f=zeros(n,1);
for i=1:m
    grad_f=grad_f+A(i,:)'/(1-A(i,:)*x);
end
for j=1:n
    grad_f(j,1)=grad_f(j,1)-(1/(1+x(j,1))-1/(1-x(j,1)));
end