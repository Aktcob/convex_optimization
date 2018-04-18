%% function
function f=func(u,W,t)
    [A B]=chol(W+diag(u));
    if B~=0
       f=inf;
    else
        f=t*ones(1,20)*u;
        f=f-log(det(W+diag(u)));
    end
end