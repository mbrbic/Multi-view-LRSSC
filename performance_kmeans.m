% Run k-means n times and report means and standard deviations of the
% performance measures.
%
% -------------------------------------------------------
% Input:
%       X:  data matrix (rows are samples)
%       k:  number of clusters
%       truth:  truth cluster indicators
%     
%
% Output:
%       CA:  clustering accuracy (mean +stdev)
%       F:  F1 measure (mean +stdev)
%       P:  precision (mean +stdev)
%       R:  recall (mean +stdev)
%       nmi:  normalized mutual information (mean +stdev)
%       AR:  adjusted rand index (mean +stdev)
%
function [CA F P R nmi AR] = performance_kmeans(X, k, truth)
max_iter = 1000; % Maximum number of iterations for KMeans
replic = 20; % Number of replications for KMeans

if (min(truth)==0)
    truth = truth+1;
end

warning('off');

for i=1:replic
    idx = kmeans(X, k,'EmptyAction','singleton','maxiter',max_iter);
    CAi(i) = 1-compute_CE(idx, truth); % clustering accuracy
    [Fi(i),Pi(i),Ri(i)] = compute_f(truth,idx); % F1, precision, recall
    nmii(i) = compute_nmi(truth,idx);
    ARi(i) = rand_index(truth,idx);
end
CA(1) = mean(CAi); CA(2) = std(CAi);
F(1) = mean(Fi); F(2) = std(Fi);
P(1) = mean(Pi); P(2) = std(Pi);
R(1) = mean(Ri); R(2) = std(Ri);
nmi(1) = mean(nmii); nmi(2) = std(nmii);
AR(1) = mean(ARi); AR(2) = std(ARi);

end
