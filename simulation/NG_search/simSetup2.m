%% Sim Setup
% addpath('C:\Program Files\MATLAB\Ipopt-3.11.8')
% addpath('C:\Program Files\MATLAB\cvx\functions\vec_') % some issue occurs when vec function is called (in my code, it happens when using log_det)
addpath('/Users/changliu/Documents/MATLAB/cvx/functions/vec_') % some issue occurs when vec function is called (in my code, it happens when using log_det)
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

sim_len = 80;
dt = 0.5;
plan_mode = 'nl'; % choose the mode of simulation: linear: use KF. nl: use gmm

tar_model = 'ped'; % static, lin, cir, sin, ped(estrian)

solver = 'sqp'; % 'sqp'

if strcmp(plan_mode,'lin')
    sensor_type = 'lin'; % rb, ran, br, lin
elseif strcmp(plan_mode,'nl')
    sensor_type = 'ran'; %rb % rb, ran, br, lin
end

inPara_sim = struct('dt',dt,'sim_len',sim_len,'sensor_type',sensor_type,'plan_mode',plan_mode);
sim = Sim(inPara_sim);

save_video = true;


%% Set field %%%
% target info
% target.pos = [17;18];%[15;15];%[30;20]; %[27;26]; %[25;35]; %[25.5;33.5]; %[25.5;30.5]; %[25.5;25.5];
target.state = [45;15;3*pi/4;1];

% linear model, used for KF
% target.A = eye(2);%[0.99 0;0 0.98];
% target.B = [0;0]; %[0.5;-0.5]; 

% nonlinear model, used for PF
switch tar_model
    case 'static'
        %%% setup for static target, KF
        target.f = @(x) x;
        target.del_f = @(x) eye(2);
%         target.A = eye(2);%[0.99 0;0 0.98];
%         target.B = [0;0]; %[0.3;-0.3];[0;0];
        target.Q = 0*eye(2); % Covariance of process noise model for the target
    
    case 'lin'
        %%% setup for moving target: linear model
        
        target.f = @(x) x+[0.5;0.5];
        target.del_f = @(x) eye(2);
        % this A, B is temporily defined to make this part compatible with KF in
        % Robot.m. Later clean this part to unify the representation of KF and PF.
        % Make sure A corresponds to del_f and B is the affine term of f.
%         target.A = eye(2);%[0.99 0;0 0.98];
%         target.B = [0.5;0.5]; %[0.5;0.5]; %[0.3;-0.3];[0;0];
        target.Q = 0.25*eye(2); %0.04 % Covariance of process noise model for the target
    
    case 'cir'
        %%% setup for moving target: circular model
        % this part is not finished yet
        des_lin_vel = 1;
        cetr = center_set(:,mode_cnt);
        center_set = [20;20]; %[[50.5;150.5],[50.5;-150.5],[-50.5;50.5],[150.5;50.5]];
        radius = norm(x-cetr);
        ang_vel = des_lin_vel/radius;
        d_ang = ang_vel*dt; % the angle increment
        cur_ang = atan2(y-cetr(2),x-cetr(1));

        tmp_x = x - radius*sin(cur_ang)*d_ang;
        tmp_y = y + radius*cos(cur_ang)*d_ang;
        %}
    
    case 'sin'
        %%% setup for moving target: sinusoidal model
        %
        u = [0.5;0.5];
        target.f = @(x) x+[u(1);u(2)*cos(0.2*x(1))];
        target.del_f = @(x) [1 0; -u(2)*sin(x(1)) 1];
%         target.A = [1 0; -u(2)*sin(x(1)) 1];
        % target.B = [0.5;0.5]; %[0.5;0.5]; %[0.3;-0.3];[0;0];
        target.Q = 0.09*eye(2); %0.04 % Covariance of process noise model for the target
        %}
    case 'ped'
        %%% setup for moving target: pedestrian (constant speed point mass) model
        % similar to robot's unicyle model except that there is no input
        % and noises are only for orientation and speed.
        target.theta_bd = [target.state(3)-5/180*pi;target.state(3)+5/180*pi];
        target.v_bd = [target.state(4)-0.5;target.state(4)+0.5];
        target.f = @(x) x+[x(4)*cos(x(3));x(4)*sin(x(3));0;0];
        target.del_f = @(x) [1 0 -x(4)*sin(x(3)) cos(x(3)); ...
            0 1 x(4)*cos(x(3)) sin(x(3));...
            0 0 1 0;...
            0 0 0 1];
        target.optf = @(x,y) x+[y(2)*cos(y(1));y(2)*sin(y(1))]; % general model f used in optimization
        target.opt_del_f = @(x,y) eye(2);
        target.Q = blkdiag(10^-1*eye(2),[0.01,0;0,0.09]);
        target.optQ = target.Q(1:2,1:2);
end

target.model_idx = 1;
target.pos_idx = [1,2]; % used for extracting the idx of pos from target state
target.traj = target.state(target.pos_idx);
target.target_model = tar_model;

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
inPara_rbt.state = [20;10;pi/3;0];%[20;20;pi/2;0];%[22;30;pi/2;0]; %[15;10;pi/2;0]; %[22;33;pi/2;0];%[40;40;pi/2;0];%;static target case:[25;15;pi/2;0];
% input constraint
inPara_rbt.a_lb = -3;
inPara_rbt.a_ub = 1;
inPara_rbt.w_lb = -pi/4;
inPara_rbt.w_ub = pi/4;
inPara_rbt.v_lb = 0;
inPara_rbt.v_ub = 3;
% robot kinematics
inPara_rbt.g = @(z,u) z+u*dt;
inPara_rbt.Qr = blkdiag(0.09*eye(2),[0.01,0;0,0.04]);
inPara_rbt.del_g = @(z,u) z+u*dt;
% target defintion
inPara_rbt.target = target;

% sensor
%%% needs further revision
inPara_rbt.sensor_type = sensor_type;
inPara_rbt.theta0 = 60/180*pi; %60/180
inPara_rbt.range = 10;%15 4.5;

if strcmp(tar_model,'ped')
    Ct = [1 0 0 0; 0 1 0 0];
else
    Ct = eye(2);
end

Cr = [1 0 0 0; 0 1 0 0];

switch sensor_type 
    case 'rb'
        % range-bearing sensor
        %     inPara_rbt.h = @(x) x-inPara_rbt.state(1:2);
        %     inPara_rbt.del_h = @(x) eye(2);
        %     inPara_rbt.R = 5*eye(2);
        inPara_rbt.h = @(x,z) x.^2-z.^2; %%%%% change this in the future. Note, R should be scaled accordingly, o.w. all particles may have very small/large weights
        inPara_rbt.del_h = @(x,z) [2*x(1) 0; 0 2*x(2)]; % z is the robot state.
        inPara_rbt.R = 400*eye(2);
        % inPara_rbt.dist_rb = 20;
    case 'ran'
        % range-only sensor
        inPara_rbt.h = @(x,z) sqrt(sum((Ct*x-z).^2)+0.1);
        inPara_rbt.opth = @(x,z) sqrt(sum((Ct(:,1:2)*x-z).^2)+0.1);
        inPara_rbt.del_h = @(x,z) Ct'*(Ct*x-z)/sqrt(sum((x-z).^2)+0.1);
        inPara_rbt.opt_del_h = @(x,z) (Ct(:,1:2)*x-z)'/sqrt(sum((Ct(:,1:2)*x-z).^2)+0.1);
        inPara_rbt.R = 1;
        inPara_rbt.mdim = 1;
    case 'br'
        % bearing-only sensor
        
    case 'lin'
        % linear sensor model for KF use
        inPara_rbt.h = @(x,z) Ct*x-z;
        inPara_rbt.opth = @(x,z) Ct(:,1:2)*x-z;
        inPara_rbt.del_h = @(x,z) Ct; % z is the robot state.
        inPara_rbt.opt_del_h = @(x,z) Ct(:,1:2);
        inPara_rbt.R = eye(2); %0.01
        inPara_rbt.mdim = 2;
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

% initializing filter
if strcmp(plan_mode,'lin')
    % KF
    inPara_rbt.est_pos = target.pos+ [5;-5];
    inPara_rbt.P = [100 0; 0 100];
    inPara_rbt.max_gmm_num = 1;
    inPara_rbt.particles = [];
elseif strcmp(plan_mode,'nl')
    % PF
    inPara_rbt.max_gmm_num = 3;
    
    if strcmp(tar_model,'ped')
        [X,Y] = meshgrid((xMin+0.5):1:(xMax-0.5),(yMin+0.5):1:(yMax-0.5));    
%         xlen = size(X,1);
%         [theta,v] = meshgrid(linspace(0,2*pi,xlen),linspace(0,2,xlen));
%         inPara_rbt.particles = [X(:),Y(:),theta(:),v(:)]';
        [theta,v] = meshgrid(linspace(target.theta_bd(1),target.theta_bd(2),5)...
            ,linspace(target.v_bd(1),target.v_bd(2),5));        
        % mesh x-y with theta-v 
        tmpxy = [X(:),Y(:)];
        tmptv = [theta(:),v(:)];
        [tmpcord1,tmpcord2] = meshgrid(1:length(tmpxy),1:length(tmptv));
        inPara_rbt.particles = [tmpxy(tmpcord1(:),:),tmptv(tmpcord2(:),:)]';
        
        inPara_rbt.est_state = target.state+[5;-5;0.3;0.5];        
        inPara_rbt.P = {}; %{[100 0; 0 100];[100 0; 0 100];[100 0; 0 100]};
    else
        [X,Y] = meshgrid((xMin+0.5):0.5:(xMax-0.5),(yMin+0.5):0.5:(yMax-0.5));
        inPara_rbt.particles = [X(:),Y(:)]'; 
        inPara_rbt.est_pos = target.pos+ [5;-5];
        inPara_rbt.P = {}; %{[100 0; 0 100];[100 0; 0 100];[100 0; 0 100]};
    end
    
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
cfg.max_gam_iter = 5; 
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
rbt = Robot7(inPara_rbt);
    
