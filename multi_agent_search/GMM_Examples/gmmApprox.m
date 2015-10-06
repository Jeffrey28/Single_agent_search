%% Chang Liu 7/3/15  changliu1289@gmail.com
%% modified for gmm approximation using weighted samples based on code online
%% algorithm and notations follow that in "EM demystified: an expectation-maximization tutorial"

%%======================================================
%% STEP 1a: Generate data from two 2D distributions.

mu1 = [1 2];       % Mean
sigma1 = [ 3 .2;   % Covariance matrix
          .2  2];
m1 = 200;          % Number of data points.

mu2 = [-1 -2];
sigma2 = [2 0;
          0 1];
m2 = 100;

% Generate sample points with the specified means and covariance matrices.
R1 = chol(sigma1);
X1 = randn(m1, 2) * R1;
X1 = X1 + repmat(mu1, size(X1, 1), 1);

R2 = chol(sigma2);
X2 = randn(m2, 2) * R2;
X2 = X2 + repmat(mu2, size(X2, 1), 1);

X = [X1; X2];

%%=====================================================
%% STEP 1b: Plot the data points and their pdfs.

figure(1);

% Display a scatter plot of the two distributions.
hold off;
plot(X1(:, 1), X1(:, 2), 'bo');
hold on;
plot(X2(:, 1), X2(:, 2), 'ro');

set(gcf,'color','white') % White background for the figure.

% First, create a [10,000 x 2] matrix 'gridX' of coordinates representing
% the input values over the grid.
gridSize = 100;
u = linspace(-6, 6, gridSize);
[A B] = meshgrid(u, u);
gridX = [A(:), B(:)];

% Calculate the Gaussian response for every value in the grid.
z1 = gaussianND(gridX, mu1, sigma1);
z2 = gaussianND(gridX, mu2, sigma2);

% Reshape the responses back into a 2D grid to be plotted with contour.
Z1 = reshape(z1, gridSize, gridSize);
Z2 = reshape(z2, gridSize, gridSize);

% Plot the contour lines to show the pdf over the data.
[C, h] = contour(u, u, Z1);
[C, h] = contour(u, u, Z2);

axis([-6 6 -6 6])
title('Original Data and PDFs');

%set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2);


%%====================================================
%% STEP 2: Choose initial values for the parameters.
% follow the initialization step in "EM demystified: an expectation-maximization tutorial"

% number of data points.
m = size(X,1);
% The vector lengths.
n = size(X,2);  
% weight for each point
v = ones(m,1);
% The number of clusters.
k = 2; 

% tolerance in the termination condition
tol = 1e-6;
%%===================================================
%% STEP 3: Run Expectation Maximization

% Loop until convergence
for iter = 1:1000
    
    fprintf('  EM Iteration %d\n', iter);

    %%===============================================
    %% STEP 3a: Expectation
    %
    if iter == 1
        % initialize gamma
        % use k-means clustering to obtain a coarse estimation of gamma
        gamma = zeros(m,n);
        idx = kmeans(X,k);
        for ii = 1:m
            gamma (ii,idx(ii)) = 1;
        end
        pre_like = -1;
    else
        % P(zi=j|yi,theta^m)
        % One row per data point, one column per cluster.
        gau_pdf = zeros(m, k);
        
        % For each cluster
        for j = 1:k
            % Evaluate the Gaussian for all data points for cluster 'j'.
            gau_pdf(:, j) = gaussianND(X, mu(j, :), sigma{j});
        end
        
        % get gamma^m_ij
        gamma = bsxfun(@times, gau_pdf, w);
        gamma = bsxfun(@rdivide, gamma, sum(gamma, 2));
        
        pre_like = sum(sum(bsxfun(@times,v,log(bsxfun(@times, gau_pdf, w)))));
    end
    %%===============================================
    %% STEP 3b: Maximization
    %%
    %% Calculate the probability for each data point for each distribution.    
    
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
    
    % Check for convergence
    for j = 1:k
        % Evaluate the Gaussian for all data points for cluster 'j'.
        gau_pdf(:, j) = gaussianND(X, mu(j, :), sigma{j});
    end
    like = sum(sum(bsxfun(@times,v,log(bsxfun(@times, gau_pdf, w)))));
    if abs(pre_like-like) < tol
        break
    end
% End of Expectation Maximization    
end

%%=====================================================
%% STEP 4: Plot the data points and their estimated pdfs.

% Display a scatter plot of the two distributions.
figure(2);
hold off;
plot(X1(:, 1), X1(:, 2), 'bo');
hold on;
plot(X2(:, 1), X2(:, 2), 'ro');

set(gcf,'color','white') % White background for the figure.

plot(mu1(1), mu1(2), 'kx');
plot(mu2(1), mu2(2), 'kx');

% First, create a [10,000 x 2] matrix 'gridX' of coordinates representing
% the input values over the grid.
gridSize = 100;
u = linspace(-6, 6, gridSize);
[A B] = meshgrid(u, u);
gridX = [A(:), B(:)];

% Calculate the Gaussian response for every value in the grid.
z1 = gaussianND(gridX, mu(1, :), sigma{1});
z2 = gaussianND(gridX, mu(2, :), sigma{2});

% Reshape the responses back into a 2D grid to be plotted with contour.
Z1 = reshape(z1, gridSize, gridSize);
Z2 = reshape(z2, gridSize, gridSize);

% Plot the contour lines to show the pdf over the data.
[C, h] = contour(u, u, Z1);
[C, h] = contour(u, u, Z2);
axis([-6 6 -6 6])

title('Original Data and Estimated PDFs');