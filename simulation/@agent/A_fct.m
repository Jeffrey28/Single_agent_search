function value = A_fct(agent,mu,sigma)
% log partition function in exponential family for Gaussian distribution
value = 1/2*mu'/sigma*mu+1/2*log(det(sigma))+log(2*pi);
end