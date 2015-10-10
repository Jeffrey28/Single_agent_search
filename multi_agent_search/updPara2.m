function [outPara] = updPara2(lambda,psi)
% This function is similar to updateProbPara. However, this one will deal
% with the case for several parameters (such as updating mu using mu_i, 
% mu_k, mu_k+1, etc...)
outPara.psi = sum(psi,3);
outPara.lambda = sum(lambda,2);
end