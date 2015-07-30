function [w,mu,sigma] = getGmmApprox(X,v,k,tol)
%% Chang Liu 7/3/15  changliu1289@gmail.com
%% gmm approximation using weighted samples
%% algorithm and notations follow the "EM demystified: an expectation-maximization tutorial"
% input: X: data points (point coordinate); v: point weight; k: number of
% clusters; tol: tolerance for termination condition
%% Choose initial values for the parameters.
% follow the initialization step in "EM demystified: an expectation-maximization tutorial"

% number of data points.
m = size(X,1);
% The vector lengths.
n = size(X,2);  

%%===================================================
%% Expectation Maximization
n_iter = 100; % maximum iteration number
% Loop until convergence
for iter = 1:n_iter    
    %%===============================================
    %% Expectation
    %
    if iter == 1
        % initialize gamma
        % use k-means clustering to obtain a coarse estimation of gamma
        gamma = zeros(m,n);
        idx = kmeans(X,k);
        for ii = 1:m
            gamma (ii,idx(ii)) = 1;
        end
        loglike = -1;
    else
        pre_loglike = loglike;
        % P(zi=j|yi,theta^m)
        % One row per data point, one column per cluster.
        gau_pdf = zeros(m, k);
        
        % For each cluster
        for j = 1:k
            % if sigma is close to singularity, make it less singular
            rcond(sigma{j})
            if (rcond(sigma{j}) < 1e-10) || (isequaln(rcond(sigma{j}),nan))
                sigma{j} = eye(2);
            end
            % Evaluate the Gaussian for all data points for cluster 'j'.
            gau_pdf(:, j) = gaussianND(X, mu(j, :), sigma{j});
        end
        
        % get gamma^m_ij
        gamma = bsxfun(@times, gau_pdf, w);
        gamma = bsxfun(@rdivide, gamma, sum(gamma, 2));
        
        loglike = sum(sum(bsxfun(@times,v,log(bsxfun(@times, gau_pdf, w)))));
        
        % check for convergence
        if abs(pre_loglike-loglike) < tol
            break
        end
    end
    %%===============================================
    %% Maximization    
    % define n^m_j = sum_i gamma^m_ij
    nj = (sum(bsxfun(@times,v,gamma),1))';
    % For each cluster
    for j = 1:k
    
        % Calculate the mode probability 
        w(j) = nj(j)/sum(v);
        
        % Calculate the mean 
        mu(j, :) = sum(bsxfun(@times,(v.*gamma(:,j)),X),1)/nj(j);
        
        % Calculate the covariance
        sigma_k = zeros(n, n);
        
        % Subtract the cluster mean from all data points.
        Xm = bsxfun(@minus, X, mu(j, :));
        
        % Calculate the contribution of each training example to the covariance matrix.
        for i = 1:m
            sigma_k = sigma_k + (v(i)*gamma(i, j)*(Xm(i, :)' * Xm(i, :)));
        end
        
        % Divide by the sum of weights.
        sigma{j} = sigma_k / nj(j);
    end
% End of Expectation Maximization    
end
% make mu a column vector
mu = mu';
fprintf('EM Iteration ends at %d\n', iter);
% this is a temporary fix for the issue of no-convergence
if iter == n_iter
    w = ones(k,1)/k;
    mu = [35;35]*ones(1,k); % assign the mean to be the target position
    for tt = 1:k
        sigma{tt} = 5*eye(2);
    end
end
end

function pdf = gaussianND(X, mu, Sigma)
% multivariate Gaussian

% Get the vector length
n = size(X, 2);

% Subtract the mean from every data point.
meanDiff = bsxfun(@minus, X, mu);

% Calculate the multivariate gaussian.
pdf = 1 / sqrt((2*pi)^n * det(Sigma)) * exp(-1/2 * sum((meanDiff/(Sigma) .* meanDiff), 2));

end

