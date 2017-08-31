% 
% Run MLRSSC on 3-sources dataset. Parameters are optimized over NMI
% measure.
%
%-------------------------------------------------------
clear;
addpath(genpath(cd))

load 3-sources
X{1} = bbc;
X{2} = guardian;
X{3} = reuters;

k = 6;
num_views = 3;
num_iter = 100;

%% Linear kernel multi-view LRSSC

fprintf('\nPairwise multiview LRSSC\n');
opts.mu = 10^2;
lambda1 = 0.3;
lambda3 = 0.3;
opts.lambda = [lambda1 (1-lambda1) lambda3];
opts.noisy = true;

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best

fprintf('\nCentroid multiview LRSSC\n');
opts.mu = 10^2;
lambda1 = 0.1;
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

opts.sigma = [50*sigma(1) 5*sigma(2) 10*sigma(3)]; 

fprintf('\nKernel pairwise multiview LRSSC\n');
opts.mu = 10^3;
lambda1 = 0.3;
lambda3 = 0.3;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best
            
fprintf('\nKernel centroid multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.3;
lambda3 = 0.7;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = centroid_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best


