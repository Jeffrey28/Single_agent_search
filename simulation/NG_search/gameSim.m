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

for ii = 1 %:sim_len
    sprintf('Progress: %d',ii/sim_len)
    %% target state update
    fld = fld.targetMove();
    
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
    
    %% robot motion planning
    %
    [optz,optu] = rbt.ngPlanner(fld);
    rbt = rbt.updState(optu);
    display('robot state:')
    display(rbt.state);
    %}
    
    % draw plot
    sim.plotFilter(rbt,fld)
%     pause()
%     sim.plotTraj(rbt,fld)
    
    
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