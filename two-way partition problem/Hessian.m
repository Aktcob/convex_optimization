%% Hessian
function H= Hessian(u,W)
    H=((W+diag(u))^-1).*((W+diag(u))^-1);
end