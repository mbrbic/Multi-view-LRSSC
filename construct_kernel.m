% Kernel construction.
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
function K = construct_kernel(X1,X2,options)

if (~exist('options','var'))
    options = [];
else
    if ~isstruct(options)
        error('parameter error!');
    end
end

if ~isfield(options,'KernelType')
    options.KernelType = 'linear';
end

switch lower(options.KernelType)
    case 'gaussian'        %  e^{-(|x-y|^2)/2t^2}
        if ~isfield(options,'sigma')
            options.sigma = 1;
        end
        if isempty(X2)
            D = eucl_dist2(X1,[],0);
        else
            D = eucl_dist2(X1,X2,0);
        end
        K = exp(-D/(2*options.sigma^2));
        
    case 'polynomial'      % (x'*y)^d
        if ~isfield(options,'d')
            options.d = 2;
        end
        if isempty(X2)
            D = full(X1 * X1');
        else
            D = full(X1 * X2');
        end
        K = D.^options.d;
        
    case 'linear'      % x'*y
        if isempty(X2)
            K = full(X1 * X1');
        else
            K = full(X1 * X2');
        end
    otherwise
        error('Parameter error: opts is not a structure.');
end

if isempty(X2)
    K = max(K,K');
end


