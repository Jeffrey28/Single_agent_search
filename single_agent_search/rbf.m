function K = rbf(u,v,l)
% rbf kernel function. 
% Input:
% u,v = matrices containing samples as rows
% l = the parameter for the kernel function
% Output:
% K = the rbf kernel matrix ( = exp(-1/(2*l^2)*(||u-v||_2)^2) )

K = zeros(size(u,1),size(v,1));
for ii = 1:size(K,1)
    for jj = 1:size(K,2)
        K(ii,jj) = exp(-1/(2*l^2)*norm(u(ii,:)-v(jj,:),2)^2);
    end
end
