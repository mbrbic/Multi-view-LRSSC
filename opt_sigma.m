function sigma = opt_sigma(X)
    N = size(X,1);
    dist = eucl_dist2(X,X); 
    dist = reshape(dist,1,N*N);
    sigma = median(dist);