%
% Solve the multiview pairwise similarities based LRSSC
%
% -------------------------------------------------------
% Input:
%       X:  cell of nxm matrices (rows are samples)
%       k:  number of clusters
%       truth:  truth cluster indicators
%       opts:  Structure value with the following fields:
%           opts.lambda:    coefficients for low-rank, sparsity and
%           consensus constraint, respectively
%           opts.num_iter:  number number of iterations
%           opts.mu:  penalty parameters for the ADMM 
%           opts.max_mu:  maximal penalty parameters for the ADMM
%           opts.rho: step size for adaptively changing mu, if 1 fixed mu
%           is used
%           opts.kernel:  kernel type (Linear, Gaussian or Polynomial)
%           opts.sigma:  parameter for Gaussian kernel
%           opts.degree:  parameter for polynomial kernel
%           opts.noisy:  true for noisy version of LRSSC false otherwise
%           opts.err_thr: convergence threshold
%
% Output:
%       Af:  joint affinity matrix
%
% version 1.0 - 20/12/2016
% version 2.0 - 30/08/2017
%
% Written by Maria Brbic (maria.brbic@irb.hr)
%
function Af = pairwise_MLRSSC(X, opts)

% setting default parameters
num_iter = 300;
mu = 100;
max_mu = 1e6;
rho = 1.5;
kernel = 'Linear';
noisy = true;
lambda = [0.3 0.3 0.3];
err_thr = 10^-3; 

if ~exist('opts', 'var')
    opts = [];
else
   if ~isstruct(opts) 
       error('Parameter error: opts is not a structure.');
   end
end

if isfield(opts, 'lambda');      lambda = opts.lambda;	end
if isfield(opts, 'num_iter');    num_iter = opts.num_iter;	end
if isfield(opts, 'mu');          mu = opts.mu;	end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;	end
if isfield(opts, 'kernel');      kernel = opts.kernel;	end
if isfield(opts, 'sigma');       sigma = opts.sigma;    end
if isfield(opts, 'degree');      degree = opts.degree;	end
if isfield(opts, 'noisy');       noisy = opts.noisy;	end
if isfield(opts, 'rho');         rho = opts.rho;	end
if isfield(opts, 'err_thr');     err_thr = opts.err_thr;	end

if strcmpi(kernel,'Gaussian') 
    noisy = true; % can't kernelize otherwise
end

mu1 = mu;
mu2 = mu;
mu3 = mu;
mu4 = mu;
mu = [mu1 mu2 mu3 mu4];

num_views = size(X,2);
n = size(X{1},1);

C1 = repmat({zeros(n,n)}, 1, num_views);
C2 = repmat({zeros(n,n)}, 1, num_views);
C3 = repmat({zeros(n,n)}, 1, num_views);
K = repmat({zeros(n,n)}, 1, num_views);
A = repmat({zeros(n,n)}, 1, num_views);
A_prev = repmat({zeros(n,n)}, 1, num_views);

Lambda1 = cell(num_views,1);
for v = 1:num_views
    Lambda1{v} = (zeros(size(X{v},2),n));
end
Lambda2 = repmat({zeros(n,n)}, 1, num_views);
Lambda3 = repmat({zeros(n,n)}, 1, num_views);
Lambda4 = repmat({zeros(n,n)}, 1, num_views);

for v = 1:num_views
    options.KernelType = kernel;
    
    if strcmpi(kernel,'Gaussian')
        options.sigma = sigma(v);
    end
    if strcmpi(kernel,'Polynomial')
        options.d = degree;
    end
    
    K{v} = construct_kernel(X{v},X{v},options);
end

iter = 0;
converged = false;

while iter < num_iter && ~converged
    
    iter = iter + 1;
    
    C_sum = repmat({zeros(n,n)}, 1, num_views);
    
    for v = 1:num_views
        for v_tmp = 1:num_views
            if v_tmp ~= v
                C_sum{v} = C_sum{v}+C2{v_tmp};
            end
        end
    end
    
    for v = 1:num_views
        A_prev{v} = A{v}; % save previous value
        [C1{v}, C2{v}, C3{v}, Lambda1{v}, Lambda2{v}, Lambda3{v}, Lambda4{v}, A{v}] = pairwise_LRSSC_1view...
            (X{v}', K{v}, num_views, C1{v}, C2{v}, C3{v}, C_sum{v}, Lambda1{v}, Lambda2{v}, Lambda3{v}, Lambda4{v}, ...
            lambda, mu, noisy);
    end
    
    % check convergence
    converged = true;
    for v=1 : num_views
        err1 = max(max(abs(A{v}-C1{v})));
        err2 = max(max(abs(A{v}-C2{v})));
        err3 = max(max(abs(A{v}-C3{v})));
        err4 = max(max(abs(A_prev{v}-A{v})));
                
        if err1>err_thr || err2>err_thr || err3>err_thr || err4>err_thr
            converged = false;
            break
        end
    end
     
    mu = min(rho*mu,max_mu); 
end

%if converged
%    fprintf('Converged in iter %d\n', iter);
%end

C_avg = zeros(n,n);
for v = 1:num_views
    C_avg = C_avg+C2{v};
end
C_avg = C_avg/(num_views);

C = C_avg;
Af = abs(C)+abs(C');

end


