function outPara = particle_filter(agent,inPara)
% particle filter given new observations
% besides, return the gmm fitting parameters.
particles = inPara.particles;
cur_pos = agent.currentPos(1:2); 
sigma_s = agent.sigma_s;
k_s = agent.k_s;
reading = inPara.obs;
n_gmm = inPara.n_gmm;
%% particle filtering
n_p = size(particles,2); % number of particles

% initalize particles weights
w = ones(n_p,1)/n_p;

% state update: since we use the static target, no update is needed

% weight update
for ii = 1:n_p
    w(ii) = sensorModel(sigma_s,k_s,particles(:,ii),cur_pos,reading);
end
w = w/sum(w);

% resampling
idx = randsample(1:n_p, n_p, true, w);
new_particles = particles(:,idx);

%% gmm fitting
tic;
% n_gmm = 10; % max cluster number
gmm_model = cell(n_gmm,1);
opt = statset('MaxIter',1000);
AIC = zeros(n_gmm,1);
for kk = 1:n_gmm
    gmm_model{kk} = fitgmdist(new_particles',kk,'Options',opt,...
        'Regularize',0.001,'CovType','full');
    AIC(kk)= gmm_model{kk}.AIC;
end
display ('gmm fitting takes time as:');
toc;

[minAIC,numComponents] = min(AIC);

best_model = gmm_model{numComponents};
gmm_mu = best_model.mu';
% gmm_sigma = zeros(2,2,numComponents);
% for ii = 1:numComponents
%     gmm_sigma(:,:,ii) = diag(best_model.Sigma(:,:,ii)); % best_model.Sigma is a vector, showing the diagonal elements of the matrix
% end
gmm_sigma = best_model.Sigma;
gmm_w = best_model.PComponents';
gmm_lambda = zeros(size(gmm_mu));
gmm_psi = zeros(size(gmm_sigma));
for jj = 1:size(gmm_mu,2)
    gmm_psi(:,:,jj) = 1/2*eye(2)/gmm_sigma(:,:,jj);
    gmm_lambda(:,jj) = gmm_sigma(:,:,jj)\gmm_mu(:,jj);
end

%% output
outPara = struct('particles',new_particles,'w',gmm_w,'mu',gmm_mu,...
    'sigma',gmm_sigma,'lambda',gmm_lambda,'psi',gmm_psi);
end

function prob = sensorModel(sigma,k,x,x_r,reading)
% x is the particle position; x_r is the robot position; reading is the
% sensor reading
pod = k/(2*pi)*(det(sigma))^(-1/2)*exp(-1/2*(x-x_r)'/sigma*(x-x_r));
if reading == 1
    prob = pod;
elseif reading == 0
    prob = 1-pod;
end
end
