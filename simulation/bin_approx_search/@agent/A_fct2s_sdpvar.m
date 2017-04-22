function value = A_fct2s_sdpvar(agent,lambda,psi)
% a short version (this is why there's a 's' in the function name)
% of the log partition function in exponential family for Gaussian distribution
% remove the log term here
% this is specially for the optimizer since a sdpvar cannot be used with
% inverse, so when psi is a sdpvar in the optimizer, I need to manually
% inverse it
inv_psi = [psi(2,2),-psi(1,2);-psi(2,1),psi(1,1)]/(psi(1,1)*psi(2,2)-psi(2,1)*psi(1,2));
value = 1/4*lambda'*inv_psi*lambda+log(pi);
end