% Euclidean Distance matrix.
%
% -------------------------------------------------------
% Input:
%       X1:  nxm data matrix (rows are samples)
%       X2:  nxm data matrix (rows are samples)
%       opts:  Structure value with the following fields:
%           opts.KernelType:    kernel type, choices are:
%           'Linear': (x'*y),'Gaussian': (e^{-(|x-y|^2)/2sigma^2}),
%           'Polynomial': ((x'*y)^d)
%           opts.sigma:  variance for Gaussian kernel
%           opts.d:  degree for polynomial kernel
%
% Output:
%       K:  nxn kernel matrix
%
%   D = EuDist(fea_a,fea_b)
%   fea_a:    nSample_a * nFeature
%   fea_b:    nSample_b * nFeature
%   D:      nSample_a * nSample_a
%       or  nSample_a * nSample_b
function [D] = eucl_dist2(X1,X2,bSqrt)

if ~exist('bSqrt','var')
    bSqrt = 1;
end

if (~exist('fea_b','var')) | isempty(X2)
    [nSmp, nFea] = size(X1);

    aa = sum(X1.*X1,2);
    ab = X1*X1';
    
    aa = full(aa);
    ab = full(ab);

    if bSqrt
        D = sqrt(repmat(aa, 1, nSmp) + repmat(aa', nSmp, 1) - 2*ab);
        D = real(D);
    else
        D = repmat(aa, 1, nSmp) + repmat(aa', nSmp, 1) - 2*ab;
    end
    
    D = max(D,D');
    D = D - diag(diag(D));
    D = abs(D);
else
    [nSmp_a, nFea] = size(X1);
    [nSmp_b, nFea] = size(X2);
    
    aa = sum(X1.*X1,2);
    bb = sum(X2.*X2,2);
    ab = X1*X2';

    aa = full(aa);
    bb = full(bb);
    ab = full(ab);

    if bSqrt
        D = sqrt(repmat(aa, 1, nSmp_b) + repmat(bb', nSmp_a, 1) - 2*ab);
        D = real(D);
    else
        D = repmat(aa, 1, nSmp_b) + repmat(bb', nSmp_a, 1) - 2*ab;
    end
    
    D = abs(D);
end

