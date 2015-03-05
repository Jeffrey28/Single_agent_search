function [outPara] = predictHumanTraj(agent,inPara)
% define input arguments
obv_traj = inPara.obv_traj; % observed trajectory of human
hor = inPara.hor;
pre_type = inPara.pre_type;
dt = inPara.mpc_dt;

% v = state([2,4]);
% initialization
pre_traj = zeros(3,hor+1); % current and predicted human future path, first row is time label
pre_traj(2:3,1) = obv_traj(2:3,end);
pre_t = (1:hor)'+obv_traj(1,end); % time points to predict at in the future
if strcmp(pre_type,'extpol')
    for ii = 1:hor
        pre_traj(:,ii+1) = pre_traj(:,ii)+v*dt;
    end
elseif strcmp(pre_type,'GP')
    obv_t = obv_traj(1,:)';
    obv_x = obv_traj(2,:)'; % observed x trajectory
    obv_y = obv_traj(3,:)';
    n = size(obv_traj,2);
    sig_n = 0; % no noise
    l = 1; % parameter for rbf
    K = rbf(obv_t,obv_t,l);
    L = chol(K+sig_n^2*eye(n),'lower');

    % predict for x-axis
    alpha_x = L'\(L\obv_x);
    k_star_x = rbf(obv_t,pre_t,l);
    x_mean = k_star_x'*alpha_x;
    v = L\k_star_x;
    x_cov = rbf(pre_t,pre_t,l)-v'*v;
    
    % predict for y-axis
    alpha_y = L'\(L\obv_y);
    k_star_y = rbf(obv_t,pre_t,l);
    y_mean = k_star_y'*alpha_y;
    v = L\k_star_y;
    y_cov = rbf(pre_t,pre_t,l)-v'*v;
end
pre_traj(:,2:end) = [pre_t;x_mean;y_mean];
outPara = pre_traj;
end