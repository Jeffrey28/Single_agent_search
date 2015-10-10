function value = A_fct2s(lambda,psi)
% a short version (this is why there's a 's' in the function name)
% of the log partition function in exponential family for Gaussian distribution
% remove the log term here
value = 1/4*lambda'/psi*lambda;
end