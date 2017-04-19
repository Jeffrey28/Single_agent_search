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
% [sim_len,sim,rbt,fld] = simSetup2();
%% %%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%% 
% optimal solution for seeding ngPlanner
optz = [];
optu = [];

for ii = 1:sim_len
    sprintf('Progress: %d',ii/sim_len)    
    
    %% target estimation
    rbt.y = rbt.sensorGen(fld);
    display('measurement')
    display(rbt.y)
%     rbt = rbt.GSF(fld);
    rbt = rbt.PF(fld);
    display('weights')
    display(rbt.wt)
    display('covariance')
    display(rbt.P);
    display('estimated position')
    display(rbt.est_pos);
    
    %% target state update
    fld = fld.targetMove();
    
    %% robot motion planning
    %
%     [optz,optu] = rbt.ngPlanner(fld,optz,optu);
    [optz,optu] = rbt.cvxPlanner(fld,optz,optu);
    rbt = rbt.updState(optu);
    display('robot state:')
    display(rbt.state);
    %}
    
    % draw plot
    sim.plotFilter(rbt,fld)
    sim.plotTraj(rbt,fld)
    pause()
    
    
    % terminating condition
%     if trace(rbt.P) <= 1 && norm(fld.target.pos-rbt.est_pos) <= 2 && norm(rbt.state(1:2)-target.pos) <= 3
%         display('target is localized')
%         break
%     end     
   %}
end

%% save simulation result
% save('test4_offset')

% run resultAnalysis.m to analyze the simulation results