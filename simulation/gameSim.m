% 9/12/16
% single-agent search. linear target motion model and linear observation
% model. 

clc 
clear % clear global variables
close all

%% Setup
% addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles');
% addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
addpath('C:\Program Files\MATLAB\Ipopt-3.11.8')
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

sim_len = 100;
dt = 1;

%%% Set Robot %%%
% Robot
inPara_rbt = struct;
inPara_rbt.state = [10;10;0;0];
inPara_rbt.sen_cov = eye(2);
inPara_rbt.inv_sen_cov = 0.02*eye(2);
inPara_rbt.max_step = sim_len;
inPara_rbt.v_lb = 0;
inPara_rbt.v_ub = 3;
inPara_rbt.w_lb = -pi/2;
inPara_rbt.w_ub = pi/2;
inPara_rbt.R = eye(2);
inPara_rbt.C = eye(2);
inPara_rbt.mpc_hor =3;
inPara_rbt.dt = dt;
inPara_rbt.theta0 = 35/180*pi;
inPara_rbt.range = 4.5;
inPara_rbt.init_pos = [50;50];
inPara_rbt.init_P = [100 0; 0 100];
inPara_rbt.gam = 1;
rbt = Robot(inPara_rbt);

%%% Set field %%%
% target info
target.pos = [30;30];
target.A = eye(2);
target.Q = eye(2); % Covariance of process noise model for the target
target.model_idx = 1;
target.traj = target.pos;

xLength = 100;%*scale; 
yLength = 100;%*scale; 
xMin = 0;
yMin = 0;
xMax = xMin+xLength;
yMax = yMin+yLength;
grid_step = 0.5; % the side length of a probability cell
inPara_fld = struct('fld_cor',[xMin,xMax,yMin,yMax],'target',target,'dt',dt);
fld = Field(inPara_fld);

% draw agents on the initial environment

%% %%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%% 
for ii = 1%:sim_len
    %% draw plot
    figure(1)
    hold on
    hdl1 = plot(rbt.traj(1,:),rbt.traj(2,:),'r','markers',5);
    set(hdl1,'MarkerFaceColor','r');
    set(hdl1,'MarkerEdgeColor','r');
    set(hdl1,'Color','r');
    set(hdl1,'LineStyle','-');
    set(hdl1,'Marker','o');    
    
    hdl2 = plot(fld.target.pos(1),fld.target.pos(2),'b','markers',5);
    set(hdl2,'MarkerFaceColor','b');
    set(hdl2,'MarkerEdgeColor','b');
    set(hdl2,'Color','b');
%     set(hdl2,'LineStyle','-');
    set(hdl2,'Marker','*');    
    
    hdl3 = plot(rbt.est_pos(1),rbt.est_pos(2),'g','markers',5);
    set(hdl3,'MarkerFaceColor','g');
    set(hdl3,'MarkerEdgeColor','g');
    set(hdl3,'Color','b');
%     set(hdl2,'LineStyle','-');
    set(hdl3,'Marker','s');    
    
    xlim([fld.fld_cor(1),fld.fld_cor(2)]);
    ylim([fld.fld_cor(3),fld.fld_cor(4)]);
    
    %% observe and update target estimation
    rbt.y = rbt.sensorGen(fld);
    rbt = rbt.KF(fld);
    
    %% robot motion planning
    [optz,optu] = rbt.Planner(fld);
    rbt = rbt.updState(optu);
end


%% save simulation result
%{
% save data
% if the data is a decimal, replace the '.' with 'p'
str_safe_dis = strrep(num2str(safe_dis),'.','p');
str_safe_marg = strrep(num2str(safe_marg),'.','p');
str_h_v = strrep(num2str(agents(1).currentV),'.','p');

str_t = datestr(clock,'dd-mmm-yyyy_HHMMSS');
folder_path = ('.\sim_res');
data_name = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.mat',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,data_name);
save(file_name,'obv_traj','est_state','pre_traj','plan_state','r_state','r_input');

% save plot
folder_path = ('.\sim_res');
fig_name = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.fig',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,fig_name);
h = gcf;
saveas (h,file_name);

% convert .fig to .pdf
fig_name2 = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.pdf',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name2 = fullfile (folder_path,fig_name2);
fig2Pdf(file_name2,300,h)
%}