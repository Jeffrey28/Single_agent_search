function idx = findEqualVec(v,M)
% find the logical array indicating if vector v is in some colunms of the matrix M
% v is n-by-1 vector and M is n-by-m matrix
len = length(v);
f = @(x,vec) (x == vec);
idx = ones(1,size(M,2));
for ii = 1:len
    idx = idx & f(v(ii),M(ii,:));
end