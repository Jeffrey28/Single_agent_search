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

% save figures to video
if save_video
    if strcmp(plan_mode,'lin')
        vidObj = VideoWriter(sprintf('search-using-KF-%s.avi',date));
    elseif strcmp(plan_mode,'nl')
        vidObj = VideoWriter(sprintf('search-using-PF-%s.avi',date));
    end
    vidObj.FrameRate = 2;
    open(vidObj);
end

for ii = 1:sim_len
    sprintf('gameSim.m, line %d, Progress: %d',MFileLineNr(),ii/sim_len)    
    
    %% target estimation
    rbt.y = rbt.sensorGen(fld);
    sprintf('gameSim.m, line %d, measurement:',MFileLineNr())
    display(rbt.y)
    
    if strcmp(plan_mode,'lin')
        rbt = rbt.KF(fld);
    elseif strcmp(plan_mode,'nl')
        %     rbt = rbt.GSF(fld);
        rbt = rbt.PF(fld);
    end
    
    sprintf('gameSim.m, line %d, weights', MFileLineNr())
    display(rbt.wt)
    sprintf('gameSim.m, line %d, covariance', MFileLineNr())
    display(rbt.P);
    sprintf('gameSim.m, line %d, estimated position', MFileLineNr())
    display(rbt.est_pos);
    
    if strcmp(plan_mode,'lin')
        sim.plotFilter_kf(rbt,fld)
    elseif strcmp(plan_mode,'nl')
        sim.plotFilter(rbt,fld)
    end
    
    %% target state update
    fld = fld.targetMove();
    
    %% robot motion planning
    %
    if strcmp(plan_mode,'lin')
        [optz,optu] = rbt.cvxPlanner_kf(fld,optz,optu);
    elseif strcmp(plan_mode,'nl')
%         [optz,optu] = rbt.ngPlanner(fld,optz,optu);
        [optz,optu] = rbt.cvxPlanner(fld,optz,optu);
    end
    
    rbt = rbt.updState(optu);
    sprintf('gameSim.m, line %d, robot state:', MFileLineNr())
    display(rbt.state);
    %}
    
    % draw plot
    sim.plotTraj(rbt,fld)
%     pause()

    % save the plot as a video
    frame_hdl = getframe(gcf);
    if save_video
        writeVideo(vidObj,frame_hdl);
    end   
    
    % terminating condition
%     if trace(rbt.P) <= 1 && norm(fld.target.pos-rbt.est_pos) <= 2 && norm(rbt.state(1:2)-target.pos) <= 3
%         display('target is localized')
%         break
%     end     
   %}
end

if save_video
    close(vidObj);
end

%% save simulation result
% save('test4_offset')

% run resultAnalysis.m to analyze the simulation results