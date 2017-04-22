function [outPara] = predictHumanTraj(agent,inPara)
% define input arguments
obv_traj = inPara.obv_traj; % observed trajectory of human
hor = inPara.hor;
pre_type = inPara.pre_type;
mpc_dt = inPara.mpc_dt;
% samp_rate = inPara.samp_rate;

% v = state([2,4]);
% initialization
pre_traj = zeros(3,hor+1); % current and predicted human future path, first row is time label
pre_traj(:,1) = obv_traj(:,end);
pre_t = (1:hor)*mpc_dt+obv_traj(1,end); % time points to predict at in the future
if strcmp(pre_type,'extpol')
    v = (obv_traj(2:3,end)-obv_traj(2:3,end-1))/(obv_traj(1,end)-obv_traj(1,end-1));
    for ii = 1:hor
        pre_traj(2:3,ii+1) = pre_traj(2:3,ii)+v*mpc_dt;
        pre_traj(1,ii+1) = pre_traj(1,ii)+mpc_dt;
    end
    outPara = struct('pre_traj',pre_traj);
elseif strcmp(pre_type,'GP')
    obv_t = obv_traj(1,:)';
    obv_x = obv_traj(2,:)'; % observed x trajectory
    obv_y = obv_traj(3,:)';
    n = size(obv_traj,2);
    sig_n = 0.01; % use a very small noise term to avoid the numerical singularity of K
    l = 1; % parameter for rbf
    K = rbf(obv_t,obv_t,l);
    L = chol(K+sig_n^2*eye(n),'lower');

    % predict for x-axis
    alpha_x = L'\(L\obv_x);
    k_star_x = rbf(obv_t,pre_t',l);
    x_mean = k_star_x'*alpha_x;
    v = L\k_star_x;
    x_cov = rbf(pre_t',pre_t',l)-v'*v;
    
    % predict for y-axis
    alpha_y = L'\(L\obv_y);
    k_star_y = rbf(obv_t,pre_t',l);
    y_mean = k_star_y'*alpha_y;
    v = L\k_star_y;
    y_cov = rbf(pre_t',pre_t',l)-v'*v;
    pre_traj(:,2:end) = [pre_t;x_mean';y_mean'];
    pre_cov = cat(3,x_cov,y_cov);
    outPara = struct('pre_traj',pre_traj,'pre_cov',pre_cov);
end
end