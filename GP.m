% A snippet for running Gaussian process regression based on the algorithm
% in Rasmussen's book
% generate data set
clear 
close all
% low = -5; high = 5;
% coe = ((high-low)*rand(5,1)+low);
% g = @(x) [x^4,x^3,x^2,x,1]*coe;
g = @(x) 3*x+1;
X = [-0.5:0.05:0.5]';
y = zeros(size(X,1),1);
for ii = 1:size(X,1)
    y(ii) = g(X(ii,:));
end
figure(1);
plot(X,y)
title('original function');
hold on
x_star = [0.6,1,2]';
y_star = g(x_star);
plot(x_star,y_star,'*','MarkerSize',6)
% get predictive mean and variance
n = size(X,1);
sig_n = 0.0001; % no noise
l = 0.1; % parameter for rbf
K = rbf(X,X,l);
L = chol(K+sig_n^2*eye(n),'lower');
alpha = L'\(L\y);
k_star = rbf(X,x_star,l);
f_mean = k_star'*alpha;
v = L\k_star;
f_cov = rbf(x_star,x_star,l)-v'*v;
plot(x_star,f_mean,'^','MarkerSize',6)