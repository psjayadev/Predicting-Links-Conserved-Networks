function [Ahat,sin_vals] = Linear_Model(X,n,N,noise_flag,Sigma_e,m)

% This function performs PCA on steady state flow measurements 

%% Getting the Constraint matrix by applying svd
if noise_flag == 0
    [u,s,~]=svd(X);
    sin_vals = s(1:n,1:n);
    Ahat = u(:,m+1:n)'; % constraint matrix
else
    L = chol(Sigma_e);
    Xs = inv(L)*X;
    [u,s,~]=svd(Xs);
    Ahat = u(:,m+1:n)'*inv(L);      % A linear model of data
    sin_vals = s(1:n,1:n);
end









