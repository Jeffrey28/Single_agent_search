%% Sim Setup
% addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles');
% addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
addpath('C:\Program Files\MATLAB\Ipopt-3.11.8')
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

sim_len = 100;
dt = 0.5;

sensor_type = 'rb'; % rb, ran, br

%% Set field %%%
% target info
target.pos = [20;30];
% linear model, used for KF
target.A = eye(2);%[0.99 0;0 0.98];
target.B = [0;0]; %[0.3;-0.3];
% nonlinear model, used for EKF
target.f = @(x) x;
target.del_f = @(x) eye(2);
target.Q = 0.01*eye(2); % Covariance of process noise model for the target
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
inPara_rbt.state = [20;20;pi/4;0];%[40;40;pi/2;0];%;
% input constraint
inPara_rbt.a_lb = -3;
inPara_rbt.a_ub = 1;
inPara_rbt.w_lb = -pi/4;
inPara_rbt.w_ub = pi/4;
inPara_rbt.v_lb = 0;
inPara_rbt.v_ub = 3;
% sensor
%%% needs further revision
inPara_rbt.sensor_type = sensor_type;
if strcmp(sensor_type,'rb')
    % range-bearing sensor
%     inPara_rbt.h = @(x) x-inPara_rbt.state(1:2);
%     inPara_rbt.del_h = @(x) eye(2);
%     inPara_rbt.R = 5*eye(2);
    inPara_rbt.h = @(x) x.^2-inPara_rbt.state(1:2).^2;
    inPara_rbt.del_h =  @(x) [2*x(1) 0; 0 2*x(2)];
    inPara_rbt.R = 5*eye(2);
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
end

inPara_rbt.theta0 = 360/180*pi;
inPara_rbt.range = 30;%4.5;

% estimation initialization
inPara_rbt.est_pos = target.pos+[5,-10,10;-0.5,10,-5];%[30;25];
inPara_rbt.P = {[100 0; 0 100];[100 0; 0 100];[100 0; 0 100]};
inPara_rbt.gmm_num = size(inPara_rbt.est_pos,2);
inPara_rbt.wt = ones(inPara_rbt.gmm_num,1)/inPara_rbt.gmm_num;

% planning
inPara_rbt.mpc_hor = 3;
inPara_rbt.dt = dt;

% simulation parameters
inPara_rbt.max_step = sim_len;
inPara_rbt.gam = 1;
rbt = Robot(inPara_rbt);