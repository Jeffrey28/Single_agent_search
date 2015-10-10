function [prob_map,gmm_data] = updateProbMap(agent,field,fit_flag)
% Construct the target distribution in the field using gaussian mixture
% mu and sigma correpsonds to the vector and matrix of mean and covariance
% of each gaussian. w is the mix
xMin = field.endpoints(1);
xMax = field.endpoints(2);
yMin = field.endpoints(3);
yMax = field.endpoints(4);
step = field.grid_step;
xLength = (xMax-xMin)/step;
yLength = (yMax-yMin)/step;
w = field.w;
mu = field.mu;
sigma = field.sigma;
% remove small weights
% max_w = max(w);
% rv_id = (abs(w) < 1e-4);
% w(rv_id) = [];
% w = w/sum(w);
% mu(:,rv_id) = [];
% sigma(:,:,rv_id) = [];

p_idx = find(w>=0); % indices of positive weight terms
n_idx = find(w<0); % indices of negative weight terms

% generate gmm model
gmm_obj_p = gmdistribution(mu(:,p_idx)',sigma(:,:,p_idx),w(p_idx)');
gmm_obj_n = gmdistribution(mu(:,n_idx)',sigma(:,:,n_idx),-w(n_idx)');

% scale factor: w in gmdistribution will be normalized. So I need to make
% scale back when summing positive and negative parts together.
scale_p = sum(w(p_idx));
scale_n = sum(-w(n_idx));

% mix_num = length(w);
[x_axis,y_axis] = meshgrid(xMin+step/2:step:xMax-step/2,yMin+step/2:step:yMax-step/2);
% prob_map = zeros(xLength,yLength);
% all_prob = zeros(numel(x_axis),mix_num);
% for ii = 1:mix_num
%     tmp = mvnpdf([x_axis(:),y_axis(:)],mu(:,ii)',sigma(:,:,ii));
%     % transform into probability mass
%     prob_map(:,:) = w(ii)*(reshape(tmp,size(prob_map))*step^2)'+prob_map(:,:);
% end

tmp = pdf(gmm_obj_p,[x_axis(:),y_axis(:)])*scale_p-pdf(gmm_obj_n,[x_axis(:),y_axis(:)])*scale_n;
prob_map = (reshape(tmp,xLength,yLength)*step^2)';

if fit_flag == 1
    % if need to generate data to fit for new gmm
    n_data = 3000;% number of randomly generated data
    gmm_data = zeros(n_data,2);
    % generate random data from gmm
    % current number of data
    num = 0;
    tic;
    while (num<n_data)
        tmp_data = random(gmm_obj_p,n_data);
        p = (scale_p*pdf(gmm_obj_p,tmp_data)-scale_n*pdf(gmm_obj_n,tmp_data))./(scale_p*pdf(gmm_obj_p,tmp_data));
        u = rand(n_data,1);
        idx = (u<=p);
        if num+sum(idx) <= n_data
            gmm_data(num+1:num+sum(idx),:) = tmp_data(idx,:);
            num = num+sum(idx);
        else
            tmp = tmp_data(idx,:);
            gmm_data(num+1:n_data,:) = tmp(1:n_data-num,:);
            num = n_data;
        end
    end
    toc;
else
    gmm_data = [];
end