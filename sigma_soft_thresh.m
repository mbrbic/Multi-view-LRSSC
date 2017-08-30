function Y=sigma_soft_thresh(X,beta)
[U S V]=svd(X,'econ');
v=soft_thresh(diag(S),beta);
Y=U*diag(v)*V';
end