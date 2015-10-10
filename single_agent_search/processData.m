% process simulation result
%% imm vs single-model
% load data
%
clear
addpath('.\sim_res');
load('h_state_1p5.mat');
sim_len = size(h_state,2); % simulation time
h_traj = h_state(1:2,:); % h_state contains [x;y;heading]
% safe_dist1 = 2; % safe_dist for imm
% safe_dist2 = 1.5; % safe_dist for extpol
% h_v = 2;

% % extrapolation
% % load('sim_traj_extpol_MPC_08-Dec-2014_1p5_2','pre_traj','r_state');
% % extpol_pre_traj-09-Dec-2014 saves the predicted trajectory using
% % extrapolation. this is a tempory way.
% load('sim_traj_extpol_greedy1_2_2_1p5_09-Dec-2014_204347','pre_traj');
% extpol_pre_traj = pre_traj(:,:,1:sim_len);
% % extpol_r_pos = r_state(1:2,1:sim_len);
% % extpol_r_v = r_state(3,1:sim_len);

% single model
load('x_pos_pre1_single_model','x_pos_pre1')
load('y_pos_pre1_single_model','y_pos_pre1')

sgm_pre_traj = zeros(2,size(x_pos_pre1,1)-1,size(x_pos_pre1,2));
sgm_pre_traj(1,:,:) = x_pos_pre1(2:end,:);
sgm_pre_traj(2,:,:) = y_pos_pre1(2:end,:);

% imm
% load('sim_traj_IMM_greedy1_2_2_1p5_09-Dec-2014_204006','pre_traj');
load('x_pos_pre_imm','x_pos_pre_imm')
load('y_pos_pre_imm','y_pos_pre_imm')

imm_pre_traj = zeros(2,size(x_pos_pre_imm,1)-1,size(x_pos_pre_imm,2));
imm_pre_traj(1,:,:) = x_pos_pre_imm(2:end,:);
imm_pre_traj(2,:,:) = y_pos_pre_imm(2:end,:);

% imm_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
% imm_r_pos = r_state(1:2,1:sim_len);
% imm_r_v = r_state(3,1:sim_len);
hor = size(imm_pre_traj,2)-1;
% hor = size(pre_traj,2)-1; % planning horizon

%% compare the prediction
pred_err = zeros(2,sim_len-1); % prediction error. 1st row for imm and 2nd row for extpol
for ii = 1:sim_len-1
    if ii+hor <= sim_len-1
        dif_vec1 = imm_pre_traj(:,:,ii)-h_traj(:,ii:ii+hor);
        dif_vec2 = sgm_pre_traj(:,:,ii)-h_traj(:,ii:ii+hor);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/(hor+1);
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/(hor+1);
%         pred_err(1,ii) = sum(sum(abs(dif_vec1),1))/(hor+1);
%         pred_err(2,ii) = sum(sum(abs(dif_vec2),1))/(hor+1);
    else
        tmp_len = sim_len-ii;
        dif_vec1 = imm_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii:sim_len-1);
        dif_vec2 = sgm_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii:sim_len-1);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/tmp_len;
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/tmp_len;
%         pred_err(1,ii) = sum(sum(abs(dif_vec1),1))/tmp_len;
%         pred_err(2,ii) = sum(sum(abs(dif_vec2),1))/tmp_len;
    end
end
ave_pred_err = mean(pred_err,2);
max_pred_err = max(pred_err,[],2);
h1 = figure;
box on
hold on
plot((1:size(pred_err,2))*0.5,pred_err(2,:),'r--','LineWidth',2)
plot((1:size(pred_err,2))*0.5,pred_err(1,:),'b','LineWidth',2)
legend('single','imm')
grid on
title('Prediction difference')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h1,'imm_vs_single','fig')
fig2Pdf('imm_vs_single',300,h1)
%}
%% compare the motion planning
% code is correct, but using the wrong data. motion planning comparison
% should compare mpc and greedy, not imm vs extrapolation
%{
pos_dif = zeros(2,sim_len-1);
v_dif = zeros(2,sim_len-1);

for ii = 1:sim_len-1
    pos_dif(1,ii) = sqrt(norm(imm_r_pos(:,ii)-h_traj(:,ii),2))-safe_dist1;
    pos_dif(2,ii) = sqrt(norm(extpol_r_pos(:,ii)-h_traj(:,ii),2))-safe_dist2;
end
ave_pos_dif = mean(pos_dif,2);
max_pos_dif = max(abs(pos_dif),[],2);
figure
hold on
plot(1:size(pos_dif,2),pos_dif(1,:),'r')
plot(1:size(pos_dif,2),pos_dif(2,:),'b')
legend('imm','extpol')
title('distance difference')

v_dif(1,:) = imm_r_v(1:end-1)-h_v;
v_dif(2,:) = extpol_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
max_v_dif = max(v_dif,[],2);
figure
hold on
plot(1:size(v_dif,2),v_dif(1,:),'r')
plot(1:size(v_dif,2),v_dif(2,:),'b')
legend('imm','extpol')
title('velocity difference')
legend('imm','extpol')
%}
%% mpc vs reactive

%% compare the motion planning
% load data
%
addpath('.\sim_res');
load('h_state_1p5.mat');
sim_len = size(h_state,2); % simulation time
h_traj = h_state(1:2,:); % h_state contains [x;y;heading]
safe_dist = 2; % safe_dist for imm
h_v = 1.5;

% mpc
load('sim_traj_IMM_MPC_2_2_1p5_17-Dec-2014_123019','pre_traj','r_state');
imm_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
imm_r_pos = r_state(1:2,1:sim_len);
imm_r_v = r_state(4,1:sim_len);
% greedy
load('sim_traj_IMM_greedy1_2_2_1p5_17-Dec-2014_132758','r_state');
extpol_pre_traj = pre_traj(:,:,1:sim_len);
extpol_r_pos = r_state(1:2,1:sim_len);
extpol_r_v = r_state(4,1:sim_len);

pos_dif = zeros(2,sim_len-1);
v_dif = zeros(2,sim_len-1);

for ii = 1:sim_len-1
    dif_vec1 = imm_r_pos(:,ii)-h_traj(:,ii);
    dif_vec2 = extpol_r_pos(:,ii)-h_traj(:,ii);
    pos_dif(1,ii) = norm(dif_vec1,2)-safe_dist;
    pos_dif(2,ii) = norm(dif_vec2,2)-safe_dist;
end
ave_pos_dif = mean(pos_dif,2);
max_pos_dif = max(abs(pos_dif),[],2);
h2 = figure;
hold on
box on
plot((1:size(pos_dif,2))*0.5,pos_dif(2,:),'r--','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(1,:),'b','LineWidth',2)
legend('reactive','mpc')
grid on
title('Distance difference')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h2,'mpc_vs_single_dis_diff','fig')
fig2Pdf('mpc_vs_single_dis_diff',300,h2)

v_dif(1,:) = imm_r_v(1:end-1)-h_v;
v_dif(2,:) = extpol_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
max_v_dif = max(abs(v_dif),[],2);
h3 = figure;
hold on
box on
plot((1:size(v_dif,2))*0.5,v_dif(2,:),'r--','LineWidth',2)
plot((1:size(v_dif,2))*0.5,v_dif(1,:),'b','LineWidth',2)
legend('reactive','mpc')
grid on
title('Velocity difference')
xlabel('time[sec]')
ylabel('speed [m/s]')
xlim([0,160])
saveas(h3,'mpc_vs_single_vel','fig')
fig2Pdf('mpc_vs_single_vel',300,h3)
 
% save('sim_res.mat','ave_pred_err','max_pred_err',...
%     'ave_pos_dif','max_pos_dif','ave_v_dif','max_v_dif');