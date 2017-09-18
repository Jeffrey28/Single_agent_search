%% Sim Setup
% addpath('C:\Program Files\MATLAB\Ipopt-3.11.8')
% addpath('C:\Program Files\MATLAB\cvx\functions\vec_') % soem issue occurs when vec function is called (in my code, it happens when using log_det)
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

sim_len = 60;
dt = 0.5;
plan_mode = 'nl'; % choose the mode of simulation: linear: use KF. nl: use gmm

if strcmp(plan_mode,'lin')
    sensor_type = 'lin'; % rb, ran, br, lin
elseif strcmp(plan_mode,'nl')
    sensor_type = 'lin'; %rb % rb, ran, br, lin
end

inPara_sim = struct('dt',dt,'sim_len',sim_len,'sensor_type',sensor_type,'plan_mode',plan_mode);
sim = Sim(inPara_sim);

save_video = true;

%% Set field %%%
% target info
target.pos = [27;26]; %[15;15]; %[25;35]; %[25.5;33.5]; %[25.5;30.5]; %[25.5;25.5];
% linear model, used for KF
target.A = eye(2);%[0.99 0;0 0.98];
target.B = [0;0]; %[0.5;-0.5]; 
% nonlinear model, used for EKF
% setup for static target, KF
target.f = @(x) x;
target.del_f = @(x) eye(2);
target.A = eye(2);%[0.99 0;0 0.98];
target.B = [0;0]; %[0.3;-0.3];[0;0];
target.Q = 0*eye(2); % Covariance of process noise model for the target

% % setup for moving target, KF
% target.f = @(x) x+[0.5;0.5];%[0.5;0.5]
% target.del_f = @(x) eye(2);
% % this A, B is temporily defined to make this part compatible with KF in
% % Robot.m. Later clean this part to unify the representation of KF and PF.
% % Make sure A corresponds to del_f and B is the affine term of f.
% target.A = eye(2);%[0.99 0;0 0.98];
% target.B = [0.5;0.5]; %[0.5;0.5]; %[0.3;-0.3];[0;0];
% target.Q = 0*eye(2); %0.04 % Covariance of process noise model for the target

target.model_idx = 1;
target.traj = target.pos;

xLength = 50;%*scale; 
yLength = 50;%*scale; 
xMin = 0;
yMin = 0;
xMax = xMin+xLength;
yMax = yMin+yLength;
grid_step = 0.5; % the side length of a probability cell
inPara_fld = struct('fld_cor',[xMin,xMax,yMin,yMax],'target',target,'dt',dt);
fld = Field(inPara_fld);

%% Set Robot %%%
% Robot
inPara_rbt = struct;
% robot state
inPara_rbt.state = [20;20;pi/2;0];%[22;30;pi/2;0]; %[15;10;pi/2;0]; %[22;33;pi/2;0];%[15;5;pi/2;0];%[40;40;pi/2;0];%;static target case:[25;15;pi/2;0];
% input constraint
inPara_rbt.a_lb = -3;
inPara_rbt.a_ub = 1;
inPara_rbt.w_lb = -pi/4;
inPara_rbt.w_ub = pi/4;
inPara_rbt.v_lb = 0;
inPara_rbt.v_ub = 3;
% robot kinematics
inPara_rbt.g = @(z,u) z+u*dt;
inPara_rbt.del_g = @(z,u) z+u*dt;
% target defintion
inPara_rbt.target = target;

% sensor
%%% needs further revision
inPara_rbt.sensor_type = sensor_type;
inPara_rbt.theta0 = 60/180*pi; %60/180
inPara_rbt.range = 10;%15 4.5;
if strcmp(sensor_type,'rb')
    % range-bearing sensor
%     inPara_rbt.h = @(x) x-inPara_rbt.state(1:2);
%     inPara_rbt.del_h = @(x) eye(2);
%     inPara_rbt.R = 5*eye(2);
    inPara_rbt.h = @(x,z) x.^2-z.^2; %%%%% change this in the future. Note, R should be scaled accordingly, o.w. all particles may have very small/large weights
    inPara_rbt.del_h = @(x,z) [2*x(1) 0; 0 2*x(2)]; % z is the robot state.
    inPara_rbt.R = 400*eye(2);
    % inPara_rbt.dist_rb = 20;
elseif strcmp(sensor_type,'ran')
    % % range-only sensor
    % inPara_rbt.C_ran = eye(2);
    % inPara_rbt.R_ran = 5;
    % inPara_rbt.dist_rb = 50;
elseif strcmp(sensor_type,'br')
    % bearing-only sensor
    % inPara_rbt.C_brg = eye(2);
    % inPara_rbt.R_brg = 0.5;
elseif strcmp(sensor_type,'lin')
    % lienar sensor model for KF use
    inPara_rbt.h = @(x,z) x-z;
    inPara_rbt.del_h = @(x,z) [1 0; 0 1]; % z is the robot state.
    inPara_rbt.R = 5*eye(2); %0.01
end

% define gamma model
% model parameters
alp1 = 1;
alp2 = 1;
alp3 = 1;
thr = 30;

inPara_rbt.alp1 = alp1; % parameters in gam
inPara_rbt.alp2 = alp2;
inPara_rbt.alp3 = alp3;
inPara_rbt.thr = thr; % the threshold to avoid very large exponential function value
inPara_rbt.tr_inc = 5; % increment factor of trust region
inPara_rbt.tr_dec = 1/2; % decrement factor of trust region
inPara_rbt.mu_inc = 3;

% estimation initialization
if strcmp(plan_mode,'lin')
    % KF
    inPara_rbt.est_pos = target.pos+ [5;-5];
    inPara_rbt.P = [100 0; 0 100];
    inPara_rbt.max_gmm_num = 1;
    inPara_rbt.particles = [];
elseif strcmp(plan_mode,'nl')    
    % xKF
    % inPara_rbt.est_pos = bsxfun(@plus,target.pos,[5,-10,10;-0.5,10,-5]);%[30;25];
    % inPara_rbt.P = {[100 0; 0 100];[100 0; 0 100];[100 0; 0 100]};
    
    % GSF
%     inPara_rbt.gmm_num = size(inPara_rbt.est_pos,2);
%     inPara_rbt.wt = ones(inPara_rbt.gmm_num,1)/inPara_rbt.gmm_num;
    % PF
    inPara_rbt.max_gmm_num = 3;%6;
    [X,Y] = meshgrid((xMin+0.5):(xMax-0.5),(yMin+0.5):(yMax-0.5));
    inPara_rbt.particles = [X(:),Y(:)]';
    inPara_rbt.est_pos = target.pos+ [5;-5];
    inPara_rbt.P = {[100 0; 0 100];[100 0; 0 100];[100 0; 0 100]};
end

% planning
inPara_rbt.mpc_hor = 3;%3;
inPara_rbt.dt = dt;

% simulation parameters
inPara_rbt.max_step = sim_len;

%%% configuration of optimization
cfg = {};

% sqp loop
cfg.initial_trust_box_size = 1;
cfg.improve_ratio_threshold = .25;
cfg.min_trust_box_size = 1e-2;%1e-2;%1e-4; % tol for sqp iteration
cfg.min_approx_improve = 1e-2;% 1e-4; % tol for sqp iteration
cfg.trust_shrink_ratio = .1; %this.tr_dec;
cfg.trust_expand_ratio = 1.5; %this.tr_inc;
cfg.max_sqp_iter = 1000; % max iter for sqp loop

% penalty loop
cfg.initial_penalty_coeff = 10;
cfg.cnt_tolerance = 1e-4; % tol for penalty iteration
cfg.max_penalty_iter = 5; %8; % max iter for penalty loop
% cfg.max_iter = 5;

% gamma loop
cfg.gamma_tol = 0.05; % tolerance for gamma iteration
cfg.max_gam_iter = 6; 
% cfg.max_merit_coeff_increases = 5;
cfg.merit_coeff_increase_ratio = 5; %10; %this.mu_inc

cfg.f_use_numerical = true;
cfg.g_use_numerical = true;
cfg.h_use_numerical = true;
cfg.full_hessian = false; %true;

cfg.del_f = target.del_f;
cfg.del_h = inPara_rbt.del_h;
cfg.callback = @(x,info) 0; % can change to plotting function later
inPara_rbt.cfg = cfg;

% inPara_rbt.gam = 1;
rbt = Robot(inPara_rbt);

