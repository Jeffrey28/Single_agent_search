function [prob_map_pf,prob_map_gmm] = updateProbMap_pf(agent,field,particles)
% Construct the target distribution in the field using particles. Also
% illustrate the map based on the fitted gmm

% mu and sigma correpsonds to the vector and matrix of mean and covariance
% of each gaussian. w is the mix
xMin = field.endpoints(1);
xMax = field.endpoints(2);
yMin = field.endpoints(3);
yMax = field.endpoints(4);
grid_step = field.grid_step;
xnum = (xMax-xMin)/grid_step; % number of points on x axis
ynum = (yMax-yMin)/grid_step;
w = field.w;
mu = field.mu;
sigma = field.sigma;

% generate the histogram-ish prob_map_pf
uni_par = (unique(particles','rows'))'; % 2-by-x matrix
n_uniq = size(uni_par,2); % # of unique points
prob_map_pf = [uni_par;zeros(3,n_uniq)]; % 4-by-n matrix, [x;y;count;clt_idx;pdf];
for ii = 1:n_uniq
    tmp = findEqualVec(uni_par(:,ii),particles);
%     tmp = abs(sum(uni_par(:,ii)*ones(1,size(particles,2))-particles,1))<1e-6;
    prob_map_pf(3,ii) = sum(tmp);
    if ~isempty(agent.clt_res)
        idx = (prob_map_pf(1,ii) == agent.clt_res(1,:))&(prob_map_pf(2,ii) == agent.clt_res(2,:));
        tmp_clt_res = agent.clt_res(3,idx); % agent.clt_res may return several '1's. 
        prob_map_pf(4,ii) = tmp_clt_res(1);
    end
end

% generate prob_map_gmm
gmm_obj = gmdistribution(mu',sigma,w');

prob_map_pf(5,:) = (pdf(gmm_obj,[prob_map_pf(1,:)',prob_map_pf(2,:)']))';

% mix_num = length(w);
[x_axis,y_axis] = meshgrid(xMin+grid_step/2:grid_step:xMax-grid_step/2,yMin+grid_step/2:grid_step:yMax-grid_step/2);

tmp = pdf(gmm_obj,[x_axis(:),y_axis(:)]);
prob_map_gmm = (reshape(tmp,xnum,ynum)*grid_step^2)';