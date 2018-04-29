clear;
%% parameters
m=300;
n=200;
A=dlmread('A.txt');
alpha=[0.01 0.01 0.2 0.2];
beta=[0.1 0.5 0.1 0.5];
for k=1:4
    disp(k);
    x=zeros(n,1);
    diff=inf;
    t=1;
    i=0;
    %% gradient method
    while(diff>=1e-3)
        i=i+1;
        f_p_diff(i)=obj_func(A,x);
        f(k,i)=f_p_diff(i);
        grad_f=grad_obj_func(A,x);
        diff_x=-grad_f;
        % find proper T
        t=1;
        while(obj_func(A,x+t*diff_x)>obj_func(A,x)+alpha(k)*t*(grad_f'*diff_x))
            t=t*beta(k);
        end
        Step_size(k,i)=t;
        % update
        x=x+t*diff_x;
        diff=norm(diff_x);
    end
end

%% plot f(x)-p
[p,I]=min(f');
f=f-p(1);
figure(1)
hold on
plot(f(1,1:I(1)),'r');
plot(f(2,1:I(2)),'b');
plot(f(3,1:I(3)),'y');
plot(f(4,1:I(4)),'g');
legend({'alpha=0.01 beta=0.1','alpha=0.01 beta=0.5','alpha=0.2 beta=0.1','alpha=0.2 beta=0.5'});
xlabel('iteration');ylabel('f-p*');
set(gca,'YScale','log');
hold off
figure(2)
plot(Step_size(2,1:I(2)),'b');xlabel('iteration');ylabel('Step size T');