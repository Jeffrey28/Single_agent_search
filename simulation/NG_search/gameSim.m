% single-agent search. linear target motion model and linear observation
% model. 
% This code is used for ACC 17
% Chang Liu 9/12/16
% The code is modified for IROS 17
% Chang Liu Jan. 2017


clc 
clear % clear global variables
close all

simSetup;
%% %%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%% 

for ii = 1:sim_len
    sprintf('Progress: %d',ii/sim_len)
    %% target state update
    fld = fld.targetMove();
    
    %% observe and update target estimation
    rbt.y = rbt.sensorGen(fld);
    display('measurement')
    display(rbt.y)
    rbt = rbt.GSF(fld);
    display('weights')
    display(rbt.wt)
    display('covariance')
    display(rbt.P);
    display('estimated position')
    display(rbt.est_pos);
    
    %% robot motion planning
    
    [optz,optu] = rbt.ngPlanner(fld);
    rbt = rbt.updState(optu);
    display('robot state:')
    display(rbt.state);
    
    % terminating condition
%     if trace(rbt.P) <= 1 && norm(fld.target.pos-rbt.est_pos) <= 2 && norm(rbt.state(1:2)-target.pos) <= 3
%         display('target is localized')
%         break
%     end     
   %}
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
% save('test4_offset')

% run resultAnalysis.m to analyze the simulation results