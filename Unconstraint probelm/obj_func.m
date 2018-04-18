function f= obj_func(A,x)
f=0;
[m,n]=size(A);
for i=1:m
    if(1-A(i,:)*x<0) 
        f=inf;
    else f=f-log(1-A(i,:)*x);
    end
end
for j=1:n
    if(1-x(j,1)*x(j,1)<0) 
        f=inf;    
    else f=f-log(1-x(j,1)*x(j,1));
    end
end