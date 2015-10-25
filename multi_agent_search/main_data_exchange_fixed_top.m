%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi Sensor-based Distributed Bayesian Estimator
% This code is for distributed bayesian estimator for target positioning and tracking
% (1) Target: Static target with unknown position
% (2) Sensors: Binary sensor with only 0 or 1 measurement
% (3) Strategy: (3.1) Observation Exchange strategy (Neighbourhood or Global-Exchange)
%               (3.2) Probability Map Consensus strategy (single step or multi-step)
%% 2015 June; All Copyright Reserved
% Fixed topology case

clear; clc; close all
Selection1 = 5;    % select observaition exchange and fusion strategies
max_EstStep = 110; % max step
switch Selection1
    case 1,  ObservExch='off'; ConsenStep=0;
    case 2,  ObservExch='off'; ConsenStep=10;
    case 3,  ObservExch='sing'; ConsenStep=0;
    case 4,  ObservExch='sing'; ConsenStep=10;
    case 5,  ObservExch='multi'; ConsenStep=0;
    case 6,  ObservExch='multi'; ConsenStep=10;
    otherwise, error('No selection.');
end

Selection2 = 1; % select the motion of agents and target
switch Selection2
    case 1,  r_move= 0; tar_move=0;
    case 2,  r_move= 0; tar_move=1;
    case 3,  r_move= 1; tar_move=0;
    case 4,  r_move= 1; tar_move=1;
    otherwise, error('No selection.');
end
% the set of robots whose pdf will be drawn
if r_move == 0
    sim_r_idx = [1,3,5];
else
    sim_r_idx = [2,4,6];
end

save_file = 0; % choose whether to save simulation results
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field Setup
fld.x = 100; fld.y = 100;  % Field size
fld.map = ones(fld.x,fld.y)/(fld.x*fld.y);
[xpt,ypt] = meshgrid(1:fld.x,1:fld.y);
fld.traj = []; % trajectory of traget

%% Target Setup
switch Selection2
    case {1,3} 
        fld.tx = 45; 
        fld.ty = 45;% Target position, but unknown to estimator
    case {2,4} 
        fld.tx = 20; 
        fld.ty = 20;% Target position, but unknown to estimator
end

fld.target.speed = 1;
fld.target.cov = 0.0001*eye(2);% covariance of the target motion model
if tar_move == 0
    fld.target.dx= 0;
    fld.target.dy= 0;
elseif tar_move == 1
    fld.target.dx= 1;
    fld.target.dy= 1;
end

%% Probability Map Consensus setup
ConsenFigure=0; % if 1, draw the concensus steps

%% Probability Map Update
% calculate the prediction matrix
% if tar_move == 1
    [ptx,pty] = meshgrid(1:fld.x,1:fld.y);
    pt = [ptx(:),pty(:)];
    upd_cell1 = cell(size(pt,1),1);
    for ii = 1:size(pt,1)
        upd_cell1{ii} = zeros(fld.x,fld.y);
        if (pt(ii,1)+fld.target.speed*fld.target.dx <= fld.x) && (pt(ii,2)+fld.target.speed*fld.target.dy <= fld.y)
            upd_cell1{ii}(pt(ii,1)+fld.target.speed*fld.target.dx,pt(ii,2)+fld.target.speed*fld.target.dy) = 1;
        end
    end
    
%% Path Planning setup
%{
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
%}

%% Multi-Robot Setup
NumOfRobot = 6;
hf1=figure(1); set(hf1,'Position',[50,50,1000,600]); % for filtering cycle
hf3 = figure(3);
if ConsenFigure==1, hCon=figure(2); set(hCon,'Position',[200,50,1000,600]); end % for consensus cycle
switch Selection2
    case {1,2,3,4}
        x_set = [20,40,60,80,60,40];
        y_set = [50,20,20,50,80,80];
end

for i=1:NumOfRobot
    if r_move == 0
        rbt(i).x = x_set(i); % sensor position.x
        rbt(i).y = y_set(i); % sensor position.x
    elseif r_move == 1
        rbt(i).T = 20; % period of circling motion
        rbt(i).center = [x_set(i);y_set(i)];
        rbt(i).r = 15;
        rbt(i).w = 2*pi/rbt(i).T;
        rbt(i).x = x_set(i);
        rbt(i).y = y_set(i)+rbt(i).r;
    end
    
    rbt(i).traj = [];
    rbt(i).map = ones(fld.x,fld.y);
    rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
%     if tar_move == 1
        rbt(i).talign_map = rbt(i).map; % store the prob_map for observations with same tiem index, i.e. P(x|z^1_1:k,...,z^N_1:k)
        rbt(i).talign_t = 0; % record the time for the time-aligned map 
%     end
    subplot(2,3,i); contourf((rbt(i).map)'); hold on; title(['Sensor ',num2str(i)]);
    rbt(i).prob = zeros(fld.x,fld.y);
    rbt(i).entropy = zeros(1,max_EstStep);
    for j = 1:NumOfRobot
        rbt(i).rbt(j).used = []; % save the observations times that have been used for updating
    end
end

% binary sensor model
sigmaVal=(fld.x/10)^2+(fld.y/10)^2; % covariance matrix for the sensor
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

%% Robot Buffer for Observation Exchange
for i=1:NumOfRobot
    for j=1:NumOfRobot
        rbtBuffer{i}.rbt(j).x=[];
        rbtBuffer{i}.rbt(j).y=[];
        rbtBuffer{i}.rbt(j).z=[];
        rbtBuffer{i}.rbt(j).k=[];
        rbtBuffer{i}.rbt(j).prob=[];
        rbtBuffer{i}.rbt(j).map = {};% record the previously calculated map to reduce computation        
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distributed Bayesian Filter for Target Tracking
count = 1; % time step
% record the time that the mpc solver cannot give a solution and the
% corresponding robot
err.time = [];
err.rbt = [];
while (1) %% Filtering Time Step
    figure(1); clf(hf1);
    
    for ii = 1:NumOfRobot
       rbt(ii).traj = [rbt(ii).traj,[rbt(ii).x;rbt(ii).y]];
    end
    fld.traj = [fld.traj,[fld.tx;fld.ty]];
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Generate measurement and observation probability
    for i=1:NumOfRobot % Robot Iteration
        rbt(i).z = sensorSim(rbt(i).x,rbt(i).y,fld.tx,fld.ty,sigmaVal);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Bayesian Updating
    % observation exchange steps:
    % (1) send/receive and update stored observations
    % (2) observe
    % (3) repeat step (1)
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
            
            %% with Multi-step observation exch
        case 'multi',
            if (Selection2 == 1) || (Selection2 == 3)
                % static target
                %% data transmission
                % (1) sending/receive   
                tempRbtBuffer = rbtBuffer;
                for i=1:NumOfRobot % Robot Iteration
                    % for information from neighbours to compare whether it is
                    % latest                   
                    for j=1:NumOfRobot
                        for t=rbt(i).neighbour
                            if t == 0 % if no neighbors
                                continue
                            end
                            % note: communication only transmit the latest
                            % observation stored in each neighbor
                            if (~isempty(rbtBuffer{t}.rbt(j).k)) && (isempty(tempRbtBuffer{i}.rbt(j).k) || (tempRbtBuffer{i}.rbt(j).k < rbtBuffer{t}.rbt(j).k))
                                tempRbtBuffer{i}.rbt(j).x = rbtBuffer{t}.rbt(j).x;
                                tempRbtBuffer{i}.rbt(j).y = rbtBuffer{t}.rbt(j).y;
                                tempRbtBuffer{i}.rbt(j).z = rbtBuffer{t}.rbt(j).z;
                                tempRbtBuffer{i}.rbt(j).k = rbtBuffer{t}.rbt(j).k;
                                tempRbtBuffer{i}.rbt(j).prob = rbtBuffer{t}.rbt(j).prob;
                            end
                        end
                    end
                end
                
                % return temperary buffer to robot buffer
                for i=1:NumOfRobot
                    for j=1:NumOfRobot
                        rbtBuffer{i}.rbt(j).x = tempRbtBuffer{i}.rbt(j).x;
                        rbtBuffer{i}.rbt(j).y = tempRbtBuffer{i}.rbt(j).y;
                        rbtBuffer{i}.rbt(j).z = tempRbtBuffer{i}.rbt(j).z;
                        rbtBuffer{i}.rbt(j).k = tempRbtBuffer{i}.rbt(j).k;
                        rbtBuffer{i}.rbt(j).prob = tempRbtBuffer{i}.rbt(j).prob;
                    end
                end
                
                % (2) observation
                % Observation of each robot
                for i=1:NumOfRobot
                        % initialize the buffer
                        rbtBuffer{i}.rbt(i).x=rbt(i).x;
                        rbtBuffer{i}.rbt(i).y=rbt(i).y;
                        rbtBuffer{i}.rbt(i).z=rbt(i).z;
                        rbtBuffer{i}.rbt(i).k=count;
                        if (~isempty(rbtBuffer{i}.rbt(i).k)) && (rbtBuffer{i}.rbt(i).z == 1)
                            rbtBuffer{i}.rbt(i).prob = sensorProb(rbtBuffer{i}.rbt(i).x,rbtBuffer{i}.rbt(i).y,fld.x,fld.y,sigmaVal);
                        elseif ~isempty(rbtBuffer{i}.rbt(i).k) && (rbtBuffer{i}.rbt(i).z == 0)
                            rbtBuffer{i}.rbt(i).prob = 1 - sensorProb(rbtBuffer{i}.rbt(i).x,rbtBuffer{i}.rbt(i).y,fld.x,fld.y,sigmaVal);
                        end
                end
                display(rbtBuffer{1}.rbt(1)),display(rbtBuffer{1}.rbt(2)),display(rbtBuffer{1}.rbt(3))
                display(rbtBuffer{1}.rbt(4)),display(rbtBuffer{1}.rbt(5)),display(rbtBuffer{1}.rbt(6))
                
                %% update by bayes rule
                % calculate probility of latest z
                for i=1:NumOfRobot % Robot Iteration
                    for j=1:NumOfRobot
                        if (~isempty(rbtBuffer{i}.rbt(j).k)) && (~ismember(rbtBuffer{i}.rbt(j).k,rbt(i).rbt(j).used))
                            rbt(i).map=rbt(i).map.*rbtBuffer{i}.rbt(j).prob;
                            rbt(i).rbt(j).used = [rbt(i).rbt(j).used,rbtBuffer{i}.rbt(j).k];
                        end
                    end
                    rbt(i).map=rbt(i).map/sum(sum(rbt(i).map));
                end               
                
            elseif (Selection2 == 2) || (Selection2 == 4)
                % moving target
                %% data transmission
                % (1) sending/receive
                % multi-step transmit of observation
                tempRbtBuffer=rbtBuffer;
                for i=1:NumOfRobot % Robot Iteration
                    % for information from neighbours to compare whether it is
                    % latest
                    for j=1:NumOfRobot
                        for t=rbt(i).neighbour
                            % note: communication only transmit the latest
                            % observation stored in each neighbor
                            if (~isempty(rbtBuffer{t}.rbt(j).k)) && (isempty(tempRbtBuffer{i}.rbt(j).k) || (tempRbtBuffer{i}.rbt(j).k(1) < rbtBuffer{t}.rbt(j).k(1)))
                                tempRbtBuffer{i}.rbt(j).x = [rbtBuffer{t}.rbt(j).x(1),tempRbtBuffer{i}.rbt(j).x];
                                tempRbtBuffer{i}.rbt(j).y = [rbtBuffer{t}.rbt(j).y(1),tempRbtBuffer{i}.rbt(j).y];
                                tempRbtBuffer{i}.rbt(j).z = [rbtBuffer{t}.rbt(j).z(1),tempRbtBuffer{i}.rbt(j).z];
                                tempRbtBuffer{i}.rbt(j).k = [rbtBuffer{t}.rbt(j).k(1),tempRbtBuffer{i}.rbt(j).k];
                            end
                        end
                        % remove unneeded records
                        rmv_idx = tempRbtBuffer{i}.rbt(j).k <= rbt(i).talign_t;
                        tempRbtBuffer{i}.rbt(j).x(rmv_idx)=[];
                        tempRbtBuffer{i}.rbt(j).y(rmv_idx)=[];
                        tempRbtBuffer{i}.rbt(j).z(rmv_idx)=[];
                        tempRbtBuffer{i}.rbt(j).k(rmv_idx)=[];
                    end
                end
                
                % return temperary buffer to robot buffer
                for i=1:NumOfRobot
                    for j=1:NumOfRobot
                        rbtBuffer{i}.rbt(j).x = tempRbtBuffer{i}.rbt(j).x;
                        rbtBuffer{i}.rbt(j).y = tempRbtBuffer{i}.rbt(j).y;
                        rbtBuffer{i}.rbt(j).z = tempRbtBuffer{i}.rbt(j).z;
                        rbtBuffer{i}.rbt(j).k = tempRbtBuffer{i}.rbt(j).k;
                    end
                end
                
                % (2) observation
                % Observation of each robot
                for i=1:NumOfRobot
                    rbtBuffer{i}.rbt(i).x=[rbt(i).x,rbtBuffer{i}.rbt(i).x];
                    rbtBuffer{i}.rbt(i).y=[rbt(i).y,rbtBuffer{i}.rbt(i).y];
                    rbtBuffer{i}.rbt(i).z=[rbt(i).z,rbtBuffer{i}.rbt(i).z];
                    rbtBuffer{i}.rbt(i).k=[count,rbtBuffer{i}.rbt(i).k];
                    % remove unneeded records
                    % only observations that are got no less than talign_t+1
                    % are used for Bayes filtering
                    rmv_idx = rbtBuffer{i}.rbt(i).k <= rbt(i).talign_t;
                    rbtBuffer{i}.rbt(i).x(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).y(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).z(rmv_idx)=[];
                    rbtBuffer{i}.rbt(i).k(rmv_idx)=[];
                end
                
                display(rbtBuffer{1}.rbt(1)),display(rbtBuffer{1}.rbt(2)),display(rbtBuffer{1}.rbt(3))
                display(rbtBuffer{1}.rbt(4)),display(rbtBuffer{1}.rbt(5)),display(rbtBuffer{1}.rbt(6))
                
                %% update by bayes rule
                
                for i=1:NumOfRobot % Robot Iteration
                    talign_flag = 1; % if all agent's observation's time are no less than talign_t+1, then talign_flag = 1, increase talign_t
                    tmp_t = rbt(i).talign_t;
                    tmp_map = rbt(i).talign_map;
                    
                    %                 if i == 1
                    %                     hf2 = figure (2);
                    %                     contourf((tmp_map)'); hold on;
                    %                     title('t_align map before prediction step');
                    %                 end
                    
                    for t = (rbt(i).talign_t+1):count
                        %% prediction step
                        tmp_map2 = zeros(size(tmp_map));
                        for k = 1:size(pt,1)
                            tmp_map2 = tmp_map2+upd_cell1{k}*tmp_map(pt(k,1),pt(k,2));
                        end
                        tmp_map = tmp_map2;
                        
                        %% updating step
                        for j=1:NumOfRobot
                            if (~isempty(rbtBuffer{i}.rbt(j).k)) && (rbtBuffer{i}.rbt(j).k(1) >= t)
                                %                             l = find(rbtBuffer{i}.rbt(j).k == t);
                                if t < rbtBuffer{i}.rbt(j).k(1)
                                    rbtBuffer{i}.rbt(j).prob = rbtBuffer{i}.rbt(j).map{t};
                                elseif t == rbtBuffer{i}.rbt(j).k(1)
                                    if rbtBuffer{i}.rbt(j).z(1) == 1
                                        rbtBuffer{i}.rbt(j).prob = sensorProb(rbtBuffer{i}.rbt(j).x(1),rbtBuffer{i}.rbt(j).y(1),fld.x,fld.y,sigmaVal);
                                    elseif rbtBuffer{i}.rbt(j).z(1) == 0
                                        rbtBuffer{i}.rbt(j).prob = 1 - sensorProb(rbtBuffer{i}.rbt(j).x(1),rbtBuffer{i}.rbt(j).y(1),fld.x,fld.y,sigmaVal);
                                    end
                                    % record the previously calculated map to
                                    % reduce computation
                                    rbtBuffer{i}.rbt(j).map{t} = rbtBuffer{i}.rbt(j).prob;
                                end
                                tmp_map = tmp_map.*rbtBuffer{i}.rbt(j).prob;
                            else
                                talign_flag = 0;
                            end
                        end
                        if (t == rbt(i).talign_t+1) && (talign_flag == 1)
                            rbt(i).talign_map = tmp_map;
                            rbt(i).talign_map = rbt(i).talign_map/sum(sum(rbt(i).talign_map));
                            tmp_t = tmp_t+1;
                        end
                    end
                    
                    rbt(i).talign_t = tmp_t;
                    rbt(i).map = tmp_map;
                    rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
                end
                % record the map for each time
                rbt(i).map_cell{count} = rbt(i).map;
            end
            %% Error
        otherwise, error('Observation Communication Error!');
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Entropy calculation of posterior pdf
    for i=1:NumOfRobot 
        tmp_map = rbt(i).map;
        % this avoids the error when some grid has zeros probability
        tmp_map(tmp_map <= realmin) = realmin;
        dis_entropy = -(tmp_map).*log2(tmp_map); % get the p*log(p) for all grid points    
%         fun = @(x,y) interp2(1:fld.x,1:fld.y,dis_entropy,x,y);
%         rbt(i).entropy(count) =  integral2(fun,1,fld.x,1,fld.y);
        rbt(i).entropy(count) = sum(sum(dis_entropy));
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot Local PDFs for Each Sensor After Including Observations
%     hf1 = figure (1); % handle for subplots 
%     for i=1:NumOfRobot 
%         subplot(2,3,i); contourf((rbt(i).map)'); hold on;
%         title(['Sensor ',num2str(i), ' Observ.= ',num2str(rbt(i).z)],'FontSize',14);
%         for j=1:NumOfRobot
%             if i==j
%                 plot(rbt(j).x, rbt(j).y, 's','Color',rbt(j).color,'MarkerSize',10,'LineWidth',3);
%             else
%                 plot(rbt(j).x, rbt(j).y, 'p','Color',rbt(j).color,'MarkerSize',10,'LineWidth',1.5);
%             end
%             plot(fld.tx, fld.ty, 'c+','MarkerSize',10,'LineWidth',3);
%             set(gca,'fontsize',14)
%         end
%     end
%     subplot(2,3,5); xlabel(['Step=',num2str(count)],'FontSize',14);
    
    % plot single figure for the first robot
    for k = sim_r_idx
        tmp_hd = figure (k+2); % handle for plot of a single robot's target PDF
        clf(tmp_hd);
        shading interp
        contourf((rbt(k).map)','LineColor','none'); 
        load('MyColorMap','mymap')
        colormap(mymap);
        colorbar
%         surfl((rbt(k).map)','light'); 
        hold on;
%         title(['Sensor ',1, ' Observ.= ',num2str(rbt(1).z)],'FontSize',16);
        for j=1:NumOfRobot
            
            % draw robot trajectory
            if j==k
                line_hdl = line(rbt(j).traj(1,:), rbt(j).traj(2,:));
                set(line_hdl,'Marker','.','Color','r','MarkerSize',3,'LineWidth',2);
                plot(rbt(j).traj(1,end), rbt(j).traj(2,end), 's','Color','r','MarkerSize',10,'LineWidth',3);
            else
                line_hdl = line(rbt(j).traj(1,:), rbt(j).traj(2,:));
                set(line_hdl,'Marker','.','Color','g','MarkerSize',3,'LineWidth',2);
                plot(rbt(j).traj(1,end), rbt(j).traj(2,end), 'p','Color','g','MarkerSize',10,'LineWidth',1.5);
            end
            
            % draw traget trajectory
            line_hdl = line(fld.traj(1,:), fld.traj(2,:));
            set(line_hdl,'Marker','.','Color','k','MarkerSize',3,'LineWidth',2);
            plot(fld.tx, fld.ty, 'k+','MarkerSize',10,'LineWidth',3);
%             plot3([fld.tx, fld.tx], [fld.ty, fld.ty],[0,0.05])
            set(gca,'fontsize',16)
        end
        xlabel(['Step=',num2str(count)],'FontSize',16);        
    end
    
    % save plots
    if (count == 1) || (count == 3) || (count == 5) || (count == 7) ||...
            (count == 10) || (count == 20) || (count == 30) || (count == 40)...
            || (count == 50) || (count == 60) || (count == 70) || (count == 80)...
            || (count == 90) || (count == 100)
        switch Selection2
            case 1,  tag = 'sta_sen_sta_tar';
            case 2,  tag = 'sta_sen_mov_tar';
            case 3,  tag = 'mov_sen_sta_tar';
            case 4,  tag = 'mov_sen_mov_tar';
        end
%         file_name1 = sprintf('./figures/data_exchange_switch/%s_%d_%s',tag,count,datestr(now,1));
%         saveas(hf1,file_name1,'fig')
%         saveas(hf1,file_name1,'jpg')
        for k = sim_r_idx
            tmp_hf = figure(k+2);
            file_name2 = sprintf('./figures/data_exchange/%s_single_%d_%d_%s',tag,k,count,datestr(now,1));
            if save_file == 1
                saveas(tmp_hf,file_name2,'fig')
                saveas(tmp_hf,file_name2,'jpg')
            end
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Probability Map Consensus
    %{
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
    %}
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% robot moves
    if r_move == 1
%         x_tmp = unidrnd(fld.x-2,1,NumOfRobot)+1; % generate random number within 2:fld.x-1
%         y_tmp = unidrnd(fld.y-2,1,NumOfRobot)+1;
        for i=1:NumOfRobot % Robot Iteration
            tmp_angl = calAngle([rbt(i).x;rbt(i).y]-rbt(i).center);
%             if ismember(i,[1,2,6])
                tmp_angl = tmp_angl+rbt(i).w;
%             elseif ismember(i,[3,4,5])
%                 tmp_angl = tmp_angl-rbt(i).w;
%             end
            rbt(i).x = rbt(i).r*cos(tmp_angl)+rbt(i).center(1);
            rbt(i).y = rbt(i).r*sin(tmp_angl)+rbt(i).center(2);
        end
    end
    
    %% Target Moves
    fld.tx = fld.tx + fld.target.speed * fld.target.dx;
    fld.ty = fld.ty + fld.target.speed * fld.target.dy;
        
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Terminate Time Cycle
    count = count+1;
    disp(count);
    if count > max_EstStep
        break
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot the entropy
hf2 = figure(2);
line_clr = ['r','g','b','c','m','k'];
line_marker = {'o','*','s','d','^','h'};
for i=1:NumOfRobot
    plot(1:count-2,rbt(i).entropy(1:count-2),line_clr(i),'LineWidth',2,'Marker',line_marker{i},'MarkerSize',2); hold on;
    xlim([0,count-1])
end
title('Entropy of the Target PDF','FontSize',16);
set(gca,'fontsize',16)
xlabel('Time','FontSize',16);
ylabel('Entropy','FontSize',16);

[hleg1, hobj1] = legend('Robot 1','Robot 2','Robot 3','Robot 4','Robot 5','Robot 6');
textobj = findobj(hobj1, 'type', 'text');
set(textobj, 'fontsize', 15);

switch Selection2
    case 1,  tag = 'sta_sen_sta_tar';
    case 2,  tag = 'sta_sen_mov_tar';
    case 3,  tag = 'mov_sen_sta_tar';
    case 4,  tag = 'mov_sen_mov_tar';
end
file_name2 = sprintf('./figures/data_exchange/%s_entropy_%s',tag,datestr(now,1));
if save_file == 1
    saveas(hf2,file_name2,'fig')
    saveas(hf2,file_name2,'jpg')
end

%% save robot structure
if save_file == 1
    file_name3 = sprintf('./figures/data_exchange/%s_robot.mat_%s',tag,datestr(now,1));
    save(file_name3,'rbt')
end