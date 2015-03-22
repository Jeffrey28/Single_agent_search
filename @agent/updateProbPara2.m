function [outPara] = updateProbPara2(agent,inPara)
% update probability map parameters based on observation results
% called before the updateProbMap.m

%% update parameters
w = inPara.w;
lambda = inPara.lambda;
psi = inPara.psi;
obs = inPara.obs;
lambda_s = inPara.lambda_s;
psi_s = inPara.psi_s;
field = inPara.field;

num = length(w);
if obs == 1 % positive observation
    mu = zeros(size(lambda));
    sigma = zeros(size(psi));
    for ii = 1:num  
        tmp_lambda = lambda_s+lambda(:,ii);
        tmp_psi = psi_s + psi(:,:,ii);
        alpha_i = A_fct2(agent,tmp_lambda,tmp_psi)-A_fct2(agent,lambda(:,ii),psi(:,:,ii))-A_fct2(agent,lambda_s,psi_s);
        w(ii) = w(ii)*exp(alpha_i);     
        lambda(:,ii) = tmp_lambda;
        psi(:,:,ii) = tmp_psi;
        sigma(:,:,ii) = 1/2*eye(2)/psi(:,:,ii);
        mu(:,ii) = sigma(:,:,ii)*lambda(:,ii);
    end
    w = w/(sum(w)); % normalize weights
elseif obs == 0 % negative observation
    w = [w;w]; % negative observation doubles the number of gaussian distribution
    lambda = [lambda,lambda];
    psi = cat(3,psi,psi);
    % update the mean and covariance for the second half parameters. 
    % The first half don't change except the weights
    k = agent.k_s;
    for ii = num+1:2*num 
        tmp_lambda = lambda_s+lambda(:,ii);
        tmp_psi = psi_s + psi(:,:,ii);
        alpha_i = A_fct2(agent,tmp_lambda,tmp_psi)-A_fct2(agent,lambda(:,ii),psi(:,:,ii))-A_fct2(agent,lambda_s,psi_s);
        
        w(ii) = -k*w(ii)*exp(alpha_i);
        lambda(:,ii) = tmp_lambda;
        psi(:,:,ii) = tmp_psi;        
    end
    w = w/sum(w);
    mu = zeros(size(lambda));
    sigma = zeros(size(psi));
    for ii = 1:2*num
        sigma(:,:,ii) = 1/2*eye(2)/psi(:,:,ii);
        mu(:,ii) = sigma(:,:,ii)*lambda(:,ii);
    end
end

%% GMM fitting
old_w = w;
old_mu = mu;
old_sigma = sigma;
old_lambda = lambda;
old_psi = psi;

tmp_field = field;
tmp_field.w = w;
tmp_field.mu = mu;
tmp_field.sigma = sigma;

% generate data and fit with a new gmm
[~,gmm_data] = agent.updateProbMap(tmp_field,1);
tic;
n_gmm = 10; % max cluster number
gmm_model = cell(n_gmm,1);
opt = statset('MaxIter',1000);
AIC = zeros(n_gmm,1);
for kk = 1:n_gmm
    gmm_model{kk} = fitgmdist(gmm_data,kk,'Options',opt,'Regularize',0.001);
    AIC(kk)= gmm_model{kk}.AIC;
end
display ('gmm fitting takes time as:');
toc;

[minAIC,numComponents] = min(AIC);

best_model = gmm_model{numComponents};
mu = best_model.mu';
sigma = best_model.Sigma;
w = best_model.PComponents';
lambda = zeros(size(mu));
psi = zeros(size(sigma));
for jj = 1:size(mu,2)
    psi(:,:,jj) = 1/2*eye(2)/sigma(:,:,jj);
    lambda(:,jj) = sigma(:,:,jj)\mu(:,jj);
end

outPara = struct('w',w,'mu',mu,'sigma',sigma,'lambda',lambda...
    ,'psi',psi,'old_w',old_w,'old_mu',old_mu,'old_sigma',old_sigma,...
    'old_lambda',old_lambda,'old_psi',old_psi);

