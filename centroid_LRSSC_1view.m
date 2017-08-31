%
% One step of ADMM for one view of centroid multi-view LRSSC.
%
% -------------------------------------------------------
% Input:
%   X: mxn data matrix for one view; columns are samples
%   K: nxn kernel matrix current view
%   C1: nxn low-rank coefficient matrix for current view
%   C2: nxn sparse coefficient matrix for current view
%   C3: nxn consensus coefficient matrix for current view
%   C_centroid: nxn centroid coefficient matrix
%   Lambda1, Lambda2, Lambda3, Lambda4:  Lagrange multipliers
%   lambda: coefficients for low-rank, sparsity and consensus constraint
%   mu: penalty parameters in augmented Lagrangian
%   noisy: if true then use noisy variant
%   opts.soft_thr:  true for soft thresholding (minimization of L1
%       norm), false for hard threhsold (minimization of L0 norm)
%
% Output:
%   C1: updated low-rank coefficient matrix for current view
%   C2: updated sparse coefficient matrix for current view
%   C3: updated consensus coefficient matrix for current view
%   Lambda1, Lambda2, Lambda3, Lambda4:  updated Lagrange multipliers
%
function [C1, C2, C3, Lambda1, Lambda2, Lambda3, Lambda4, A] = centroid_LRSSC_1view(X, K, C1, C2, C3, C_centroid, ...
    Lambda1, Lambda2, Lambda3, Lambda4, lambda, mu, noisy)

n = size(X,2);

% updating A
if ~noisy
    inv_A = inv(mu(1)*K+(mu(2)+mu(3)+mu(4))*eye(n));
    A = inv_A * (mu(1)*K+mu(2)*(C2-Lambda2/mu(2))+mu(3)*(C1-Lambda3/mu(3))+mu(4)*(C3-Lambda4/mu(4))+X'*Lambda1);
else
    inv_A = inv(K+(mu(2)+mu(3)+mu(4))*eye(n));
    A = inv_A * (K+mu(2)*(C2-Lambda2/mu(2))+mu(3)*(C1-Lambda3/mu(3))+mu(4)*(C3-Lambda4/mu(4)));
end

% updating C1,C2 and C3
C2_new = soft_thresh(A+Lambda2/mu(2),lambda(2)/mu(2));
C2_new = C2_new - diag(diag(C2_new));
C1_new = sigma_soft_thresh(A+Lambda3/mu(3),lambda(1)/mu(3));
C3_new = 1/(2*lambda(3)+mu(4))*(2*lambda(3)*C_centroid+mu(4)*A+Lambda4);

%update variables
C1 = C1_new;
C2 = C2_new;
C3 = C3_new;

% updating Lagrange multipliers
if ~noisy
    Lambda1 = Lambda1 + mu(1) * (X - X*A);
end
Lambda2 = Lambda2 + mu(2) * (A - C2);
Lambda3 = Lambda3 + mu(3) * (A - C1);
Lambda4 = Lambda4 + mu(4) * (A - C3);

end

