%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi Sensor-based Distributed Bayesian Estimator
% This code is for distributed bayesian estimator for target positioning and tracking
% (1) Target: Static target with unknown position
% (2) Sensors: Binary sensor with only 0 or 1 measurement
% (3) Strategy: (3.1) Observation Exchange strategy (Neighbourhood or Global-Exchange)
%               (3.2) Probability Map Consensus strategy (single step or multi-step) 
%% 2015 June; All Copyright Reserved

clear; clc; close all
Selection = 5;    % select strategies
max_EstStep = 100; % max step
switch Selection
    case 1,  ObservExch='off'; ConsenStep=0; 
    case 2,  ObservExch='off'; ConsenStep=10;
    case 3,  ObservExch='sing'; ConsenStep=0;
    case 4,  ObservExch='sing'; ConsenStep=10;
    case 5,  ObservExch='multi'; ConsenStep=0;
    case 6,  ObservExch='multi'; ConsenStep=10;
    otherwise, error('No selection.');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Field Setup
fld.x = 50; fld.y = 50;  % Field size
fld.map = ones(fld.x,fld.y)/(fld.x*fld.y);

%% Target Steup
fld.tx = 35; fld.ty = 20;% Target position, but unknown to estimator
fld.target.speed=0; 
fld.target.dx=-0.5;
fld.target.dy=0.7;

%% Binary Sensor Model
sigmaVal=(fld.x/3)^2+(fld.y/3)^2;

%% Probability Map Consensus setup
ConsenFigure=0;
% ConsenStep=0; 

%% Observation Exchange strategy setup
% ObservExch='multi';  % 'off', 'sing', 'multi'

%% Multi-Robot Setup
NumOfRobot = 6; 
hFig=figure(1); set(hFig,'Position',[50,50,1000,600]); % for filtering cycle
if ConsenFigure==1, hCon=figure(2); set(hCon,'Position',[200,50,1000,600]); end % for consensus cycle
for i=1:NumOfRobot
    rbt(i).x = 1 + fld.x/(1+NumOfRobot)*i; % sensor position.x
    rbt(i).y = 1 ; % sensor position.x
    rbt(i).speed = 1;
    rbt(i).map = ones(fld.x,fld.y); 
    rbt(i).map = rbt(i).map/sum(sum(rbt(i).map));
    subplot(2,3,i); contourf(rbt(i).map); hold on; title(['Sensor ',num2str(i)]);
    rbt(i).prob = zeros(fld.x,fld.y);
end
%% Communication structure
rbt(1).neighbour=[2,6];
rbt(2).neighbour=[1,3];
rbt(3).neighbour=[2,4];
rbt(4).neighbour=[3,5];
rbt(5).neighbour=[4,6];
rbt(6).neighbour=[5,1];

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
while (1) %% Filtering Time Step
    figure(1);clf(hFig);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Target Movement
    fld.tx = fld.tx + fld.target.speed * fld.target.dx; 
    fld.ty = fld.ty + fld.target.speed * fld.target.dy; 
    
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
            rbtBuffer{i}.rbt(i).x=rbt(i).x;
            rbtBuffer{i}.rbt(i).y=rbt(i).y;
            rbtBuffer{i}.rbt(i).z=rbt(i).z;
            rbtBuffer{i}.rbt(i).k=count; 
        end
        % multi-step transmit of observation
        for i=1:NumOfRobot % Robot Iteration
            % for information from neighbours to compare whether it is
            % latest
            tempRbtBuffer{i}=rbtBuffer{i};
            for j=1:NumOfRobot
                    for t=rbt(i).neighbour                        
                    if (tempRbtBuffer{i}.rbt(j).k <= rbtBuffer{t}.rbt(j).k)
                        tempRbtBuffer{i}.rbt(j).x = rbtBuffer{t}.rbt(j).x;
                        tempRbtBuffer{i}.rbt(j).y = rbtBuffer{t}.rbt(j).y;
                        tempRbtBuffer{i}.rbt(j).z = rbtBuffer{t}.rbt(j).z;
                        tempRbtBuffer{i}.rbt(j).k = rbtBuffer{t}.rbt(j).k;
                    end
                end
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
        
        % update by bayse rule
        for i=1:NumOfRobot % Robot Iteration
            for j=1:NumOfRobot      
                rbt(i).map=rbt(i).map.*rbtBuffer{i}.rbt(j).prob;                               
            end
            rbt(i).map=rbt(i).map/sum(sum(rbt(i).map));
        end
        %% Error
        otherwise, error('Observation Communication Error!');     
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot Local PDFs for Each Sensor
    for i=1:NumOfRobot
        subplot(2,3,i); contourf(rbt(i).map); hold on;
        title(['Sensor ',num2str(i), ' Observ.= ',num2str(rbt(i).z)]);
        for j=1:NumOfRobot
        if i==j
            plot(rbt(j).x, rbt(j).y, 'rs','MarkerSize',8,'LineWidth',3);
        else
            plot(rbt(j).x, rbt(j).y, 'mp','MarkerSize',8,'LineWidth',1.5);
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
        for i=1:NumOfRobot
            rbtCon(i).map=tempRbtCon(i).map;
            if ConsenFigure==1
                subplot(2,3,i); contourf(rbtCon(i).map); title(['Sensor ',num2str(i)]);
                hold on;
                for j=1:NumOfRobot
                if i==j
                    plot(rbt(j).x, rbt(j).y, 'rs', 'MarkerSize',8,'LineWidth',3);
                else
                    plot(rbt(j).x, rbt(j).y, 'mp', 'MarkerSize',8,'LineWidth',1.5);
                end
            end
            end
        end
    end
    for i=1:NumOfRobot % Robot Iteration
        rbt(i).map=rbtCon(i).map;
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Path Planning for Multi-Robots
    for i=1:NumOfRobot
        % find peak PDF
        maxValue=max(max(rbt(i).map));
        [colOfPeak,  rowOfPeak] = find(rbt(i).map == maxValue)
        % calculate moving direction in unity
        direction.x=rowOfPeak - rbt(i).x;
        direction.y=colOfPeak - rbt(i).y;
        lengthXY=sqrt(direction.x^2 + direction.y^2);
        direction.x=direction.x/lengthXY;
        direction.y=direction.y/lengthXY;
        % move to peak PDF
%         rbt(i).x = rbt(i).x + rbt(i).speed * direction.x; 
%         rbt(i).y = rbt(i).y + rbt(i).speed * direction.y;
        DistCheck=0;
        for t=rbt(i).neighbour
            DistWithNeigh = sqrt((rbt(i).x - rbt(t).x)^2 + (rbt(i).y - rbt(t).y)^2);
            if DistWithNeigh < 5 % Distance with neighbours are too short
               DistCheck = DistCheck + 1;  
            end
        end        
        if (DistCheck  <= 0) % Enough distance
            rbt(i).x = rbt(i).x + rbt(i).speed * direction.x; 
            rbt(i).y = rbt(i).y + rbt(i).speed * direction.y;
        else % Too short distance--> go to random direction
            rbt(i).x = rbt(i).x - rbt(i).speed * (rand(1)-0.5); 
            rbt(i).y = rbt(i).y - rbt(i).speed * (rand(1)-0.5);
        end        
     end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Terminate Time Cycle
    count = count+1;
    disp(count);
    if count > max_EstStep
        break
    end
end