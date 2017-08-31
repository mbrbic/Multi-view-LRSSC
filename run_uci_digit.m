% 
% Run MLRSSC on UCI digit dataset. Parameters are optimized over NMI
% measure.
%
%-------------------------------------------------------
clear;
addpath(genpath(cd))

load uci-digit
X{1} = mfeat_fac;
X{2} = mfeat_fou;
X{3} = mfeat_kar;

k = 10;
num_views = 3;
num_iter = 100;

%% Linear kernel multi-view LRSSC

fprintf('\nPairwise multiview LRSSC\n');
opts.mu = 10^2;
lambda1 = 0.5;
lambda3 = 0.7;
opts.lambda = [lambda1 (1-lambda1) lambda3];
opts.noisy = false;

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best

fprintf('\nCentroid multiview LRSSC\n');
opts.mu = 10;
lambda1 = 0.9;
lambda3 = 0.7;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = centroid_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best

%% Gaussian kernel multi-view LRSSC

opts.kernel = 'Gaussian';
opts.err_thr = 10^(-5);
for v=1:num_views
   sigma(v) = opt_sigma(X{v});
end
opts.sigma = [5*sigma(1) 0.5*sigma(2) 0.5*sigma(3)]; 

fprintf('\nKernel pairwise multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.7;
lambda3 = 0.7;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best
            
fprintf('\nKernel centroid multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.7;
lambda3 = 0.5;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = centroid_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best

