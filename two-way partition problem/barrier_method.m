%% barrier method
clear;
%% initial parameters
alpha=0.01;
beta=0.5;
W=dlmread('W.txt');
Size=size(W);
[V,D]=eig(W);
u(:,1)=(1-D(1))*ones(Size(1),1);
t(1)=1;
i=1;
while(1)
    while(1)
        f=func(u(:,i),W,t);
        D_f(i,1)=-f/t;
        grad_f=grad_func(u(:,i),W,t);
        H=Hessian(u(:,i),W);
        diff_u=-H\grad_f;
        lamda=grad_f'*(H^-1)*grad_f;
        if(lamda<=2*1e-8) 
            break;
        end
        % find proper T
        T=1;
        while(func(u(:,i)+T*diff_u,W,t)>func(u(:,i),W,t)+alpha*T*(grad_f'*diff_u))
            T=T*beta;
        end
        % update
        u(:,i)=u(:,i)+T*diff_u;
    end
    u(:,i+1)=u(:,i);
    i=i+1;
    t=t*15;
    if (Size(1)/t<=1e-6)
        break;
    end
end
d_star=-sum(u(:,i))
D_f=D_f-d_star;
plot(D_f);;xlabel('iteration');ylabel('-f/t-p*');
% set(gca,'YScale','log');
X_star=15*(W+diag(u(:,i)))^-1/t;
rank_X=rank(X_star)