%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi Sensor-based Distributed Bayesian Estimator
% This code is for distributed bayesian estimator for target positioning and tracking
% (1) Target: Static target with unknown position
% (2) Sensors: Binary sensor with only 0 or 1 measurement
% (3) Strategy: (3.1) Observation Exchange strategy (Neighbourhood or Global-Exchange)
%               (3.2) Probability Map Consensus strategy (single step or multi-step)
%% 2015 June; All Copyright Reserved

clear; clc; close all
Selection1 = 5;    % select observaition exchange and fusion strategies
max_EstStep = 100; % max step
switch Selection1
    case 1,  ObservExch='off'; ConsenStep=0;
    case 2,  ObservExch='off'; ConsenStep=10;
    case 3,  ObservExch='sing'; ConsenStep=0;
    case 4,  ObservExch='sing'; ConsenStep=10;
    case 5,  ObservExch='multi'; ConsenStep=0;
    case 6,  ObservExch='multi'; ConsenStep=10;
    otherwise, error('No selection.');
end

Selection2 = 2; % select the motion of agents and target
switch Selection2
    case 1,  r_move= 0; tar_move=0;
    case 2,  r_move= 0; tar_move=1;
    case 3,  r_move= 1; tar_move=0;
    case 4,  r_move= 1; tar_move=1;
    otherwise, error('No selection.');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field Setup
fld.x = 50; fld.y = 50;  % Field size
fld.map = ones(fld.x,fld.y)/(fld.x*fld.y);
[xpt,ypt] = meshgrid(1:fld.x,1:fld.y);

%% Target Steup
fld.tx = 15; fld.ty = 15;% Target position, but unknown to estimator
fld.target.speed = 1;
fld.target.cov = 5*eye(2);% covariance of the target motion model
if tar_move == 0
    fld.target.dx= 0;
    fld.target.dy= 0;
elseif tar_move == 1
    fld.target.dx= 0.5;
    fld.target.dy= 0.5;
end

%% Probability Map Consensus setup
ConsenFigure=0; % if 1, draw the concensus steps

%% Probability Map Update
% calculate the prediction matrix
if tar_move == 1
    [ptx,pty] = meshgrid(1:fld.x,1:fld.y);
    pt = [ptx(:),pty(:)];
    upd_cell1 = cell(size(pt,1),1);
    for ii = 1:size(pt,1)
        upd_cell1{ii} = mvnpdf(pt,pt(ii,:)+[fld.target.speed*fld.target.dx,fld.target.speed*fld.target.dy],fld.target.cov);
        upd_cell1{ii} = reshape(upd_cell1{ii},fld.x,fld.y);
    end
    
    upd_cell2 = cell(size(pt,1),1);
    for ii = 1:size(pt,1)
        upd_cell2{ii} = mvnpdf(pt,pt(ii,:)-[fld.target.speed*fld.target.dx,fld.target.speed*fld.target.dy],fld.target.cov);
        upd_cell2{ii} = reshape(upd_cell2{ii},fld.x,fld.y);
    end
end

%% Path Planning setup
hor = 3; % planning horizon
% desired distances among neighboring agents
% desDist = 10*[0 1 sqrt(3) 0 sqrt(3) 1; 
%     1 0 1 sqrt(3) 0 sqrt(3); 
%     sqrt(3) 1 0 1 sqrt(3) 0; 
%     0 sqrt(3) 1 0 1 sqrt(3); 
%     sqrt(3) 0 sqrt(3) 1 0 1; 
%     1 sqrt(3) 0 sqrt(3) 1 0];
desDist = 10*[0 1 sqrt(3) 2 sqrt(3) 1; 
    1 0 1 sqrt(3) 2 sqrt(3); 
    sqrt(3) 1 0 1 sqrt(3) 2; 
    2 sqrt(3) 1 0 1 sqrt(3); 
    sqrt(3) 2 sqrt(3) 1 0 1; 
    1 sqrt(3) 2 sqrt(3) 1 0];

dt = 1; % discretization time
vl = 0; vu = 3; % lower and upper bounds for robot speed

usePathUpdMap = 0; % if 1, each robot update its map based on neighbors' planned path; if 0, no update is conducted
useSimObs = 0; % if 1, each robot simulates the observations for neighbors' planned path

% precompute combinatorial matrix for calculating POD
all_comb = {};
for ii = 1:hor
    all_comb = [all_comb;num2cell(nchoosek(1:hor,ii),2)];
end
%% Multi-Robot Setup
NumOfRobot = 6;
hf1=figure(1); set(hf1,'Position',[50,50,1000,600]); % for filtering cycle
if ConsenFigure==1, hCon=figure(2); set(hCon,'Position',[200,50,1000,600]); end % for consensus cycle
% x_set = 10*([0,sqrt(3)/2,sqrt(3)/2,0,-sqrt(3)/2,-sqrt(3)/2]+1);
% y_set = 10*([1,1/2,-1/2,-1,-1/2,1/2]+2);
x_set = [5:5:30];
y_set = 10*ones(1,NumOfRobot);
for i=1:NumOfRobot
    rbt(i).x = x_set(i)+0.1*rand(1,1); % sensor position.x
    rbt(i).y = y_set(i)+0.1*rand(1,1); % sensor position.x
    rbt(i).speedLimit = 1;
    rbt(i).map = ones(fld.x,fld.y);
    rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
    if tar_move == 1
        rbt(i).talign_map = rbt(i).map;
        rbt(i).talign_t = 0;
    end
    subplot(2,3,i); contourf((rbt(i).map)'); hold on; title(['Sensor ',num2str(i)]);
    rbt(i).prob = zeros(fld.x,fld.y);
end

% binary sensor model
sigmaVal=(fld.x/5)^2+(fld.y/5)^2; % covariance matrix for the sensor
k_s = 2*pi*sqrt(det(sigmaVal)); % normalizing factor
s_psi = 1/2*eye(2)/sigmaVal; % psi for the sensor
% robot colors
rbt(1).color = 'r';
rbt(2).color = 'g';
rbt(3).color = 'y';
rbt(4).color = 'c';
rbt(5).color = 'm';
rbt(6).color = 'w';

%% Communication structure
rbt(1).neighbour=[2,6];
rbt(2).neighbour=[1,3];
rbt(3).neighbour=[2,4];
rbt(4).neighbour=[3,5];
rbt(5).neighbour=[4,6];
rbt(6).neighbour=[1,5];
% rbt(1).neighbour=[2,3,4,5,6];
% rbt(2).neighbour=[1,3,4,5,6];
% rbt(3).neighbour=[1,2,4,5,6];
% rbt(4).neighbour=[1,2,3,5,6];
% rbt(5).neighbour=[1,2,3,4,6];
% rbt(6).neighbour=[1,2,3,4,5];
%% Robot Buffer for Observation Exchange
for i=1:NumOfRobot
    for j=1:NumOfRobot
        rbtBuffer{i}.rbt(j).x=rbt(j).x;
        rbtBuffer{i}.rbt(j).y=rbt(j).y;
        rbtBuffer{i}.rbt(j).z=0;
        rbtBuffer{i}.rbt(j).k=0;
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distributed Bayesian Filter for Target Tracking
count = 1;
% record the time that the mpc solver cannot give a solution and the
% corresponding robot
err.time = [];
err.rbt = [];
while (1) %% Filtering Time Step
    figure(1); clf(hf1);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Target Movement
    if rem(count,2) == 1
        fld.tx = fld.tx + fld.target.speed * fld.target.dx;
        fld.ty = fld.ty + fld.target.speed * fld.target.dy;
    else
        fld.tx = fld.tx - fld.target.speed * fld.target.dx;
        fld.ty = fld.ty - fld.target.speed * fld.target.dy;
    end
    
    %% Generate measurement and observation probability
    for i=1:NumOfRobot % Robot Iteration
        rbt(i).z = sensorSim(rbt(i).x,rbt(i).y,fld.tx,fld.ty,sigmaVal);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bayesian Updating
    switch ObservExch
        %% No observation exch
        case 'off',
            for i=1:NumOfRobot % Robot Iteration
                if rbt(i).z == 1
                    rbt(i).prob = sensorProb(rbt(i).x,rbt(i).y,fld.x,fld.y,sigmaVal);
                elseif rbt(i).z == 0
                    rbt(i).prob = 1 - sensorProb(rbt(i).x,rbt(i).y,fld.x,fld.y,sigmaVal);
                end
                rbt(i).map = rbt(i).map.*rbt(i).prob;
                rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
            end
            
            %% With single observation exch
        case 'sing',
            for i=1:NumOfRobot % Robot Iteration
                if rbt(i).z == 1
                    rbt(i).prob = sensorProb(rbt(i).x,rbt(i).y,fld.x,fld.y,sigmaVal);
                elseif rbt(i).z == 0
                    rbt(i).prob = 1 - sensorProb(rbt(i).x,rbt(i).y,fld.x,fld.y,sigmaVal);
                end
            end
            for i=1:NumOfRobot
                rbt(i).map = rbt(i).map.*rbt(i).prob;
                for t=rbt(i).neighbour
                    rbt(i).map = rbt(i).map.*rbt(t).prob;
                end
                rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
            end
            
            %% with Multis-step observation exch
        case 'multi',
            % Observation of each robot
            for i=1:NumOfRobot
                if count == 1
                    % initialize the buffer
                    rbtBuffer{i}.rbt(i).x=rbt(i).x;
                    rbtBuffer{i}.rbt(i).y=rbt(i).y;
                    rbtBuffer{i}.rbt(i).z=rbt(i).z;
                    rbtBuffer{i}.rbt(i).k=count;
                elseif count > 1
                    rbtBuffer{i}.rbt(i).x=[rbt(i).x,rbtBuffer{i}.rbt(i).x];
                    rbtBuffer{i}.rbt(i).y=[rbt(i).y,rbtBuffer{i}.rbt(i).y];
                    rbtBuffer{i}.rbt(i).z=[rbt(i).z,rbtBuffer{i}.rbt(i).z];
                    rbtBuffer{i}.rbt(i).k=[count,rbtBuffer{i}.rbt(i).k];
                    % remove unneeded records
                    rmv_idx = rbtBuffer{i}.rbt(i).k < rbt(i).talign_t;
                    rbtBuffer{i}.rbt(i).x(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).y(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).z(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).k(rmv_idx)=[];
                end
            end
            
            % multi-step transmit of observation
            for i=1:NumOfRobot % Robot Iteration
                % for information from neighbours to compare whether it is
                % latest
                tempRbtBuffer{i}=rbtBuffer{i};
                for j=1:NumOfRobot
                    for t=rbt(i).neighbour
                        % note: communication only transmit the latest
                        % observation stored in each neighbor
                        if (tempRbtBuffer{i}.rbt(j).k(1) <= rbtBuffer{t}.rbt(j).k(1))
                            tempRbtBuffer{i}.rbt(j).x = [rbtBuffer{t}.rbt(j).x(1),tempRbtBuffer{i}.rbt(j).x];
                            tempRbtBuffer{i}.rbt(j).y = [rbtBuffer{t}.rbt(j).y(1),tempRbtBuffer{i}.rbt(j).y];
                            tempRbtBuffer{i}.rbt(j).z = [rbtBuffer{t}.rbt(j).z(1),tempRbtBuffer{i}.rbt(j).z];
                            tempRbtBuffer{i}.rbt(j).k = [rbtBuffer{t}.rbt(j).k(1),tempRbtBuffer{i}.rbt(j).k];
                        end
                    end
                     % remove unneeded records
                     rmv_idx = tempRbtBuffer{i}.rbt(j).k < rbt(i).talign_t;
                     tempRbtBuffer{i}.rbt(j).x(rmv_idx)=[];
                     tempRbtBuffer{i}.rbt(j).y(rmv_idx)=[];
                     tempRbtBuffer{i}.rbt(j).z(rmv_idx)=[];
                     tempRbtBuffer{i}.rbt(j).k(rmv_idx)=[];
                end
            end
            % return temperary buffer to robut buffer
            for i=1:NumOfRobot
                for j=1:NumOfRobot
                    rbtBuffer{i}.rbt(j).x = tempRbtBuffer{i}.rbt(j).x;
                    rbtBuffer{i}.rbt(j).y = tempRbtBuffer{i}.rbt(j).y;
                    rbtBuffer{i}.rbt(j).z = tempRbtBuffer{i}.rbt(j).z;
                    rbtBuffer{i}.rbt(j).k = tempRbtBuffer{i}.rbt(j).k;
                end
            end
            
            % calculate probility of latest z
            for i=1:NumOfRobot % Robot Iteration
                for j=1:NumOfRobot
                    if rbtBuffer{i}.rbt(j).z == 1
                        rbtBuffer{i}.rbt(j).prob = sensorProb(rbtBuffer{i}.rbt(j).x,rbtBuffer{i}.rbt(j).y,fld.x,fld.y,sigmaVal);
                    elseif rbtBuffer{i}.rbt(j).z == 0
                        rbtBuffer{i}.rbt(j).prob = 1 - sensorProb(rbtBuffer{i}.rbt(j).x,rbtBuffer{i}.rbt(j).y,fld.x,fld.y,sigmaVal);
                    end
                end
            end
            
            % update by bayes rule
            
            for i=1:NumOfRobot % Robot Iteration
                % prediciton step
                if tar_move == 1
                    [ptx,pty] = meshgrid(1:fld.x,1:fld.y);
                    pt = [ptx(:),pty(:)];
                    tmp_map = zeros(size(rbt(i).map));
                    
                    
                    if rem(count,2) == 1
                        for k = 1:size(pt,1)
                            tmp_map = tmp_map+upd_cell1{k}*rbt(i).map(pt(k,1),pt(k,2));
                        end
                    else
                        for k = 1:size(pt,1)
                            tmp_map = tmp_map+upd_cell2{k}*rbt(i).map(pt(k,1),pt(k,2));
                        end
                    end
                    rbt(i).map = tmp_map;
                end
                
                % updating step
                for j=1:NumOfRobot
                    rbt(i).map=rbt(i).map.*rbtBuffer{i}.rbt(j).prob;
                end
                rbt(i).map=rbt(i).map/sum(sum(rbt(i).map));
            end
            %% Error
        otherwise, error('Observation Communication Error!');
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot Local PDFs for Each Sensor After Including Observations
    hf1 = figure (1); 
%     clf(hf1)
    for i=1:NumOfRobot 
        subplot(2,3,i); contourf((rbt(i).map)'); hold on;
        title(['Sensor ',num2str(i), ' Observ.= ',num2str(rbt(i).z)]);
        for j=1:NumOfRobot
            if i==j
                plot(rbt(j).x, rbt(j).y, 's','Color',rbt(j).color,'MarkerSize',8,'LineWidth',3);
            else
                plot(rbt(j).x, rbt(j).y, 'p','Color',rbt(j).color,'MarkerSize',8,'LineWidth',1.5);
            end
            plot(fld.tx, fld.ty, 'c+','MarkerSize',8,'LineWidth',3);
        end
    end
    subplot(2,3,5); xlabel(['Step=',num2str(count)]);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Probability Map Consensus
    for i=1:NumOfRobot % Robot Iteration
        rbtCon(i).map=rbt(i).map;
    end
    for ConStep=1:ConsenStep % Consensus cycle
        if ConsenFigure==1, figure(2);clf(hCon); end
        for i=1:NumOfRobot % Robot Iteration
            neighNum=length(rbt(i).neighbour)+1;
            for t=rbt(i).neighbour
                tempRbtCon(i).map=rbtCon(i).map+rbtCon(t).map;
            end
            tempRbtCon(i).map=tempRbtCon(i).map/neighNum;
        end
        % plot local PDFs after concensus
        for i=1:NumOfRobot
            rbtCon(i).map=tempRbtCon(i).map;
            if ConsenFigure==1
                subplot(2,3,i); contourf((rbtCon(i).map)'); title(['Sensor ',num2str(i)]);
                hold on;
                for j=1:NumOfRobot
                    if i==j
                        plot(rbt(j).x, rbt(j).y, 's','Color',rbt(j).color,'MarkerSize',8,'LineWidth',3);
                    else
                        plot(rbt(j).x, rbt(j).y, 'p','Color',rbt(j).color, 'MarkerSize',8,'LineWidth',1.5);
                    end
                end
            end
        end
    end
    for i=1:NumOfRobot % Robot Iteration
        rbt(i).map=rbtCon(i).map;
    end
         
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% robot moves
    if r_move == 1
        x_tmp = unidrnd(fld.x,1,NumOfRobot);
        y_tmp = unidrnd(fld.y,1,NumOfRobot);
        for i=1:NumOfRobot % Robot Iteration
            rbt(i).x = x_tmp(i);
            rbt(i).y = y_tmp(i);
        end
    end
    
    % save plots
    %{
    file_name1 = sprintf('./figures/5connect/sim_0.3_0.1_1_1_actual_peak/fig1_%d',count);
    saveas(hf1,file_name1,'jpg')
    file_name2 = sprintf('./figures/5connect/sim_0.3_0.1_1_1_actual_peak/fig3_%d',count);
    saveas(hf3,file_name2,'jpg')
    %}
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Terminate Time Cycle
    count = count+1;
    disp(count);
    if count > max_EstStep
        break
    end
    
end