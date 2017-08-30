% 
% Run MLRSSC on Reuters digit dataset. Parameters are optimal over NMI
% measure.
%
%-------------------------------------------------------
clear;
addpath(genpath(cd))

num_views = 5;
k = 6;

X{1} = spconvert(load('Index_EN-EN_sample'));
X{2} = spconvert(load('Index_EN-FR_sample'));
X{3} = spconvert(load('Index_EN-GR_sample'));
X{4} = spconvert(load('Index_EN-IT_sample'));
X{5} = spconvert(load('Index_EN-SP_sample'));

truth = load('reuters_truth');

num_iter = 100;

%% Linear kernel multi-view LRSSC

fprintf('\nPairwise multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.9;
lambda3 = 0.7;
opts.lambda = [lambda1 (1-lambda1) lambda3];
opts.noisy = true;

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best

fprintf('\nCentroid multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.5;
lambda3 = 0.3;
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

opts.sigma = [50*sigma(1) 5*sigma(2) 10*sigma(3)]; % FIXME

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


