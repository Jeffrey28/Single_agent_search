% single-agent search. linear target motion model and linear observation
% model. 
% This code is used for ACC 17
% Chang Liu 9/12/16
% The code is modified for IROS 17 (not submitted)
% Chang Liu Jan. 2017
% The code continues being modified for general NGP
% Chang Liu Apr. 2017

% clc 
clear % clear global variables
close all

simSetup;
dbstop if error
% [sim_len,sim,rbt,fld] = simSetup2();
%% %%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%% 
% optimal solution for seeding ngPlanner
optz = [];
optu = [];

merit_set = zeros(3,sim_len); % used for debugging the merit function values. debug purpose only

% save figures to video
%%% in the future, write a function that can detect duplicate file name and
%%% automatically assign a new name
if save_video
    if strcmp(plan_mode,'lin')
        vidObj = VideoWriter(sprintf('search-using-KF-static-%s.avi',date));
    elseif strcmp(plan_mode,'nl')
        vidObj = VideoWriter(sprintf('search-using-PF-%s.avi',date));
    end
    vidObj.FrameRate = 2;
    open(vidObj);
end

tic
for ii = 1:sim_len
    fprintf('[main loop] gameSim.m, line %d, iteration %d, Progress: %d\n',MFileLineNr(),ii,ii/sim_len)    
    
    %% target moves
    fld = fld.targetMove();
    
    %% target estimation
    rbt.y = rbt.sensorGen(fld);
%     sprintf('gameSim.m, line %d, measurement:',MFileLineNr())
%     display(rbt.y)
    
    if strcmp(plan_mode,'lin')
        rbt = rbt.KF(fld);
    elseif strcmp(plan_mode,'nl')
        %     rbt = rbt.GSF(fld);
        rbt = rbt.PF(fld);
    end
    
%     sprintf('gameSim.m, line %d, weights', MFileLineNr())
%     display(rbt.wt)
%     sprintf('gameSim.m, line %d, covariance', MFileLineNr())
%     display(rbt.P);
%     sprintf('gameSim.m, line %d, estimated position', MFileLineNr())
%     display(rbt.est_pos);
    
    if strcmp(plan_mode,'lin')
        sim.plotFilter_kf(rbt,fld)
    elseif strcmp(plan_mode,'nl')
        sim.plotFilter(rbt,fld)
    end
   
    %% robot motion planning
    %
    if strcmp(plan_mode,'lin')
%         [optz,optu] = rbt.cvxPlanner_kf(fld,optz,optu);
%         [optz,optu,s,snum,merit, model_merit, new_merit] = rbt.cvxPlanner_scp(fld,optz,optu,plan_mode);
        [optz,optu,s,snum,merit, model_merit, new_merit] = rbt.cvxPlanner_ipopt(fld,optz,optu,plan_mode);
    elseif strcmp(plan_mode,'nl')
%         [optz,optu] = rbt.ngPlanner(fld,optz,optu);
%         [optz,optu] = rbt.cvxPlanner(fld,optz,optu);
%         [optz,optu] = rbt.cvxPlanner_sqp(fld,optz,optu);
        [optz,optu,s,snum,merit, model_merit, new_merit] = rbt.cvxPlanner_scp(fld,optz,optu,plan_mode);
%         [optz,optu,s,snum,merit, model_merit, new_merit] = rbt.cvxPlanner_ipopt(fld,optz,optu,plan_mode);
    end
    
    merit_set(:,ii) = [merit;model_merit;new_merit];
    
    rbt = rbt.updState(optu);
    rbt.snum = snum;
%     fprintf('[main loop] gameSim.m, line %d, robot state:\n', MFileLineNr())
%     display(rbt.state);
    %
    
    % debug purpose only
    % show the solution of path planner. Note that z in this solution is
    % different from optz since optz is computed using kinematic model as
    % optu (which is same as u here).
    
    u = rbt.convState(s,snum,'u');
    z = rbt.convState(s,snum,'z');
    x = rbt.convState(s,snum,'x'); 
    P = rbt.convState(s,snum,'P');
    x_pred = rbt.convState(s,snum,'x_pred');
    P_pred = rbt.convState(s,snum,'P_pred');
    K = rbt.convState(s,snum,'K');
    %}
    
    % draw plot
    sim.plotTraj(rbt,fld)
%     pause(0.5)

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
fprintf('gameSim.m, Line %d. gmm fitting takes time as:',MFileLineNr());
toc

if save_video
    close(vidObj);
end

%% save simulation result
% save('test4_offset')

% run resultAnalysis.m to analyze the simulation results