function value = A_fct2(agent,lambda,psi)
% log partition function in exponential family for Gaussian distribution
value = 1/4*lambda'/psi*lambda-1/2*log(det(psi))+log(pi);
end