function [w,mu,sigma] = updateProbPara(w,mu,sigma,obs,mu_s,sigma_s)
num = length(w);
if obs == 1 % positive observation
    for ii = 1:num
        alpha_i = 1;
        w(ii) = w(ii)*exp(alpha_i);
        tmp_sigma = inv(sigma_s)+inv(sigma(:,:,ii));
        tmp_mu = tmp_sigma*(mu_s/sigma_s+mu(ii)/sigma(:,:,ii));
        mu(ii) = tmp_mu;
        sigma(ii) = tmp_sigma;
    end
    w = w/(sum(w)); % normalize weights
elseif obs == 0 % negative observation
    w = [w;w]; % negative observation doubles the number of gaussian distribution
    % update the mean and covariance for the second half parameters. 
    % The first half don't change except the weights
    for ii = num+1:2*num 
        tmp_sigma = inv(sigma_s)+inv(sigma(:,:,ii));
        tmp_mu = tmp_sigma*(mu_s/sigma_s+mu(ii)/sigma(:,:,ii));
        mu(ii) = tmp_mu;
        sigma(ii) = tmp_sigma;
        alpha_i = 1;
        w(ii) = -w(ii)*alpha_i;
    end
    w = w/sum(w);
end
end

function A = logPartition(mu,sigma)
% may be able to simplify if we use diagonal covariance matrix for sensor
% model
1/2*mu*sigma\mu
end    