clear;
%% parameters
m=300;
n=200;
alpha=0.01;
beta=0.5;
A=dlmread('A.txt');
x=zeros(n,1);
f=obj_func(A,x);
lamda=inf;
t=1;
i=0;
%% gradient method
while(1)
    i=i+1;
    f(i)=obj_func(A,x);
    grad_f=grad_obj_func(A,x);
    H=Hessian(A,x);
    diff_x=-H\grad_f;
    lamda=grad_f'*(H^-1)*grad_f;
    if(lamda<=2*1e-8) 
        break;
    end
    % find proper T
    t=1;
    while(obj_func(A,x+t*diff_x)>obj_func(A,x)+alpha*t*(grad_f'*diff_x))
        t=t*beta;
    end
    T(i)=t;
    % update
    x=x+t*diff_x;
    disp(lamda);
end
%% plot f(x)-p
[~,P]=size(f);
p=f(P)
f=f-p;
figure(1)
plot(f,'r');xlabel('iteration');ylabel('f-p*');
figure(2)
plot(T,'b');xlabel('iteration');ylabel('Step size');