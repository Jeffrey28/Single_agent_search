function prob = calGmmProb(w,mu,sig,x)
%% calculate the probability of given points x, using the gmm
prob = 0;
for ii = 1:length(w)
    prob = prob+w(ii)*exp(-(x-mu(:,ii))'/sig{ii}*(x-mu(:,ii))/2)/(2*pi*sqrt(det(sig{ii})));
end