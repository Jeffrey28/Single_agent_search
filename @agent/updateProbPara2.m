function [w,mu,sigma,lambda,psi] = updateProbPara2(agent,w,lambda,psi,obs,lambda_s,psi_s)
% update probability map parameters based on observation results
% called before the updateProbMap.m
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
end