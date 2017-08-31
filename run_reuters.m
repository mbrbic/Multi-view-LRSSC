% 
% Run MLRSSC on Reuters digit dataset. Parameters are optimized over NMI
% measure.
%
%-------------------------------------------------------
clear;
addpath(genpath(cd))

load reuters

X{1} = spconvert(EN_EN_sample);
X{2} = spconvert(EN_FR_sample);
X{3} = spconvert(EN_GR_sample);
X{4} = spconvert(EN_IT_sample);
X{5} = spconvert(EN_SP_sample);

k = 6;
num_iter = 100;
num_views = 5;

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

opts.sigma = [5*sigma(1) 50*sigma(2) 50*sigma(3) 50*sigma(4) 10*sigma(5)];

fprintf('\nKernel pairwise multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.5;
lambda3 = 0.9;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = pairwise_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best
            
fprintf('\nKernel centroid multiview LRSSC\n');
opts.mu = 10^4;
lambda1 = 0.5;
lambda3 = 0.3;
opts.lambda = [lambda1 (1-lambda1) lambda3];

A = centroid_MLRSSC(X, opts); % joint affinity matrix
[best.CA best.F best.P best.R best.nmi best.AR] = spectral_clustering(A, k, truth);
best


