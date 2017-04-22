% single-agent search. linear target motion model and linear observation
% model. 
% This code is used for ACC 17
% Chang Liu 9/12/16

clc 
clear % clear global variables
close all

%% Setup
% addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles');
% addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
addpath('C:\Program Files\MATLAB\Ipopt-3.11.8')
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

sim_len = 50;
dt = 0.5;

%%% Set field %%%
% target info
target.pos = [20;30];
target.A = eye(2);%[0.99 0;0 0.98];
target.B = [0.3;-0.3];
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

%%% Set Robot %%%
% Robot
inPara_rbt = struct;
inPara_rbt.state = [10;10;pi/4;0];%[40;40;pi/2;0];%;
inPara_rbt.sen_cov = eye(2);
inPara_rbt.max_step = sim_len;
inPara_rbt.a_lb = -3;
inPara_rbt.a_ub = 1;
inPara_rbt.w_lb = -pi/4;
inPara_rbt.w_ub = pi/4;
inPara_rbt.v_lb = 0;
inPara_rbt.v_ub = 3;
inPara_rbt.R = eye(2);
inPara_rbt.C = eye(2);
inPara_rbt.mpc_hor = 3;
inPara_rbt.dt = dt;
inPara_rbt.theta0 = 60/180*pi;
inPara_rbt.range = 5;%4.5;
inPara_rbt.est_pos = target.pos+[0.5;-0.5];%[30;25];
inPara_rbt.P = [10 0; 0 10];
inPara_rbt.gam = 1;
rbt = Robot(inPara_rbt);

%% %%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%% 
for ii = 1:sim_len
    sprintf('Progress: %d',ii/sim_len)
    %% target state update
    fld = fld.targetMove();
    
    %% observe and update target estimation
    rbt.y = rbt.sensorGen(fld);
    display('measurement')
    display(rbt.y)
    rbt = rbt.KF(fld);
    display('covariance')
    display(rbt.P);
    display('estimated position')
    display(rbt.est_pos);
    
    %% robot motion planning
    [optz,optu] = rbt.Planner2(fld)
    rbt = rbt.updState(optu);
    display('robot state:')
    display(rbt.state);
    
    if trace(rbt.P) <= 1 && norm(fld.target.pos-rbt.est_pos) <= 2 && norm(rbt.state(1:2)-target.pos) <= 3
        display('target is localized')
        break
    end        
end

%% draw plot
figure
hold on
hdl1 = plot(rbt.traj(1,:),rbt.traj(2,:),'r','markers',3);
set(hdl1,'MarkerFaceColor','r');
set(hdl1,'MarkerEdgeColor','r');
set(hdl1,'Color','r');
set(hdl1,'LineStyle','-');
set(hdl1,'Marker','o');

hdl2 = plot(fld.target.traj(1,:),fld.target.traj(2,:),'b','markers',3);
set(hdl2,'MarkerFaceColor','b');
set(hdl2,'MarkerEdgeColor','b');%
set(hdl2,'Color','b');
%     set(hdl2,'LineStyle','-');
set(hdl2,'Marker','*');

% draw FOV
a1 = rbt.traj(3,end)-rbt.theta0;  % A random direction
a2 = rbt.traj(3,end)+rbt.theta0;
t = linspace(a1,a2,50);
x0 = rbt.traj(1,end);
y0 = rbt.traj(2,end);
x1 = rbt.traj(1,end) + rbt.r*cos(t);
y1 = rbt.traj(2,end) + rbt.r*sin(t);
plot([x0,x1,x0],[y0,y1,y0],'y-','LineWidth',1.5)

legend('robot','target');%,'estimated target')

xlim([fld.fld_cor(1),fld.fld_cor(2)]);
ylim([fld.fld_cor(3),fld.fld_cor(4)]);
box on
axis equal

%% save simulation result
save('test4_offset')

% run resultAnalysis.m to analyze the simulation results