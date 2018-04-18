%% grad_func
function grad_f=grad_func(u,W,t)
    Size=size(W);
    grad_f=t*ones(Size(1),1)-diag((W+diag(u))^-1);
end