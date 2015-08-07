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
[xpt,ypt] = meshgrid(1:fld.x,1:fld.y);

%% Target Steup
fld.tx = 35; fld.ty = 35;% Target position, but unknown to estimator
fld.target.speed=0;
fld.target.dx=-0.5;
fld.target.dy=0.7;

%% Binary Sensor Model


%% Probability Map Consensus setup
ConsenFigure=0; % if 1, draw the concensus steps
% ConsenStep=0;

%% Observation Exchange strategy setup
% ObservExch='multi';  % 'off', 'sing', 'multi'

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
% rbt(1).neighbour=[2,3,5,6];
% rbt(2).neighbour=[1,3,4,6];
% rbt(3).neighbour=[1,2,4,5];
% rbt(4).neighbour=[2,3,5,6];
% rbt(5).neighbour=[1,3,4,6];
% rbt(6).neighbour=[1,2,4,5];
rbt(1).neighbour=[2,3,4,5,6];
rbt(2).neighbour=[1,3,4,5,6];
rbt(3).neighbour=[1,2,4,5,6];
rbt(4).neighbour=[1,2,3,5,6];
rbt(5).neighbour=[1,2,3,4,6];
rbt(6).neighbour=[1,2,3,4,5];
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
     
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Path Planning for Multi-Robots
    % transmit planned path to neighbors
    if count == 1
        % initialize the plan_path
        for i=1:NumOfRobot
            for j=1:NumOfRobot
                rbtBuffer{i}.rbt(j).plan_path = [rbt(j).x;rbt(j).y]*ones(1,hor+1);% the last one is predicted by the robot motion model
            end
        end
    end
    
    if count ~= 1
        % transmit
        for i=1:NumOfRobot
            for j=1:NumOfRobot
                if ismember(j,rbt(i).neighbour)
                    rbtBuffer{i}.rbt(j).plan_path = rbt(j).plan_path;
                end
            end
        end
    end
    
    % define the map used for path planning: plan_map
    for i=1:NumOfRobot      
        rbt(i).plan_map = rbt(i).map;
    end
    % update plan_map based on neighbors's planned path
    if usePathUpdMap == 1 && count ~= 1
        for i=1:NumOfRobot
            for t=rbt(i).neighbour
                tmp_path = rbtBuffer{i}.rbt(t).plan_path;
                tmp_prob = 1;
                for kk = 1:size(tmp_path,2)
                    if useSimObs == 1
                        % simulate the observation result p(z=1)=\sum
                        % p(z|x^t)p(x^t)
                        % probability of detection
                        d_prob = sum(sum(rbt(i).plan_map.*sensorProb(tmp_path(1,kk),tmp_path(2,kk),fld.x,fld.y,sigmaVal)));
                        % simulate observation in planning horizon
                        tmp_obs = (rand(1,1) < d_prob);
                        
                        if tmp_obs == 0
                            tmp_prob = tmp_prob*(1 - sensorProb(tmp_path(1,kk),tmp_path(2,kk),fld.x,fld.y,sigmaVal));
                        elseif tmp_obs == 1
                            tmp_prob = tmp_prob*sensorProb(tmp_path(1,kk),tmp_path(2,kk),fld.x,fld.y,sigmaVal);
                        end
                    elseif useSimObs == 0
                        tmp_prob = tmp_prob*(1 - sensorProb(tmp_path(1,kk),tmp_path(2,kk),fld.x,fld.y,sigmaVal));
                    end
                end
            end
            rbt(i).plan_map = rbt(i).plan_map*tmp_prob;
            rbt(i).plan_map = rbt(i).plan_map/sum(sum(rbt(i).plan_map));
        end
    end
    
    % path planning for each robot    
    for i=1:NumOfRobot 
        % find peak PDF
        maxValue=max(max(rbt(i).plan_map));
        [xPeak, yPeak] = find(rbt(i).plan_map == maxValue);
%         xPeak = fld.tx;
%         yPeak = fld.ty;
        rbt(i).peak(:,count) = [xPeak(1);yPeak(1)];
        % define varaibles for Yalmip
        x = sdpvar(2,hor+1); % robot state, the current state is the first column [x;y]
%         u = sdpvar(2,hor); % robot input; [v;theta]
        u = sdpvar(2,hor); % robot input; [dx;dy]
        
        %% define objective function
        % find initial solution
        %{
        init_x = zeros(size(x));
        init_u = zeros(size(u));
        init_x(:,1) = [rbt(i).x;rbt(i).y];
        for kk = 1:hor
            ang = calAngle([xPeak(1);yPeak(1)]-init_x(:,kk));
            init_u(kk) = ang;
            init_x(:,kk+1) = init_x(:,kk)+rbt(i).speed*[cos(init_u(kk));sin(init_u(kk))];
        end
%         assign(x,init_x);
%         assign(u,init_u);
        %}
        
        %% obj1: distance to peak position (choose one if there are multiple)
%         obj1 = sum((x(1:2,2)-[xPeak(1);yPeak(1)]).^2);
        obj1 = 0;
        
        for t = rbt(i).neighbour
           tmp_plan_path = rbtBuffer{i}.rbt(t).plan_path;
           tmp_vpeak = x(:,2:end) - rbt(i).peak(:,count)*ones(1,hor);
           for ll = 1:hor    
               obj1 = obj1+(tmp_vpeak(1,ll)^2+tmp_vpeak(2,ll)^2)^2;
           end
        end
%         obj1 = 0;
        
        %% obj2: probability mass
        obj2 = 0;
        %
        % fit the grid using gmm
        % need to check how matlab decides the k
        % another issue: it's possible for the fitting to diverge. May need
        % to add certain minimum gmm_w,mu,sig to deal with this problem. Or
        % fix this issue by improving this fitting method.
        k = 1;
        tol = 1e-5;
        [gmm_w,gmm_mu,gmm_sig] = getGmmApprox([xpt(:),ypt(:)],reshape(rbt(i).plan_map',numel(rbt(i).plan_map),1),k,tol);
        
        
        % objective function to guide the robot towards high probability
        % position at every step
        %{
        for k = 1:hor
            obj2 = obj2+1-calGmmProb(gmm_w,gmm_mu,gmm_sig,x(1:2,k+1));
        end
        %}
        
        % objective function to guide the robot to follow the probability
        % of high info gain
        % A funcitons [A_1,...,A_hor], A_1,...A_hor are for future
        % positions, i.e. A(Phi^R_(k+ii)) in (27b) in DSCC2015 paper
        s_Af = sdpvar(hor,1); 
        for ii = 1:hor
            s_Af(ii) = A_fct2s(sigmaVal\x(1:2,ii+1),s_psi);
        end
        
        % lambda and psi for the target
        gmm_lambda = zeros(size(gmm_mu));
        gmm_psi = cell(size(gmm_sig));
        for j = 1:length(gmm_w)
            alpha = sdpvar(size(all_comb,1),1);
            gmm_lambda(:,j) = gmm_sig{j}\gmm_mu(:,j);
            gmm_psi{j} = 1/2*eye(2)/gmm_sig{j};
            gmm_Af = A_fct2s(gmm_lambda(:,j),gmm_psi{j});
            
            tmp_obj = 0;
            tmp_psi_prod = gmm_psi{j}; % this term is used for representing the constant term in exp(alpha)
            for kk = 1:size(all_comb,1) % give the term in the form of exp(A_(ij))-exp(A_i)-exp(A_j)+1
                tmp_idx = all_comb{kk};
                
                % assemble all lambda's for calculating sum_lambda -- the sum of
                % all lambda's
                s_lambda = sigmaVal\x(1:2,tmp_idx+1); % get the corresponding lambda for the robot
                tmp_lambda_set = [gmm_lambda(:,j),s_lambda];

                % assemble all psi's for calculating tmp_psi -- the sum of
                % all psi's
                tmp_psi_set = zeros(2,2,length(tmp_idx)+1);
                tmp_psi_set(:,:,1) = gmm_psi{j};
                for ll = 1:length(tmp_idx)
                    tmp_psi_set(:,:,ll+1) = s_psi;
                end
                
                sum_para = updPara2(tmp_lambda_set,tmp_psi_set);
                sum_psi = sum_para.psi;
                sum_lambda = sum_para.lambda;
                %         for ll = 1:length(tmp_idx)
                tmp_psi_prod = tmp_psi_prod*(2*s_psi)^length(tmp_idx);
                %         end
                tmp_psi_prod = tmp_psi_prod/sum_psi;
                % it seems that if the robot is stuck in a small region for a few
                % steps, the tmp_psi_prod tends to be singular. Haven't looked into
                % why this happens
%                 if abs(det(tmp_psi_prod)) < 1e-6
%                     tmp_psi_prod = tmp_psi_prod+eye(2)*1e-6;
%                 end
                tmp_coef = sqrt(det(tmp_psi_prod))*(2*pi)^(-length(tmp_idx));
                alpha(kk) = A_fct2s(sum_lambda,sum_psi)-gmm_Af-sum(s_Af(tmp_idx));
%                 alpha(kk) = 1;
%                 tmp_obj = tmp_obj+(-1)^length(tmp_idx)*(alpha(kk));
                tmp_obj = tmp_obj+(-1)^length(tmp_idx)*k_s^length(tmp_idx)*tmp_coef*exp(alpha(kk));
            end
            
            %
            if is(tmp_obj,'complex')
                display('complex obj')
            end
            %}
            
            obj2 = obj2+gmm_w(j)*(1+tmp_obj);
        end
        %}
        
        %% obj3: distance among neighbors
%         obj3 = 0;
        %
%         peak = [fld.tx;fld.ty];
        tmp_obj = 0;
        for t = rbt(i).neighbour
           tmp_plan_path = rbtBuffer{i}.rbt(t).plan_path;
           tmp_v = x(:,2:end) - tmp_plan_path(:,2:end);
%            tmp_vpeak = x(:,2:end) - peak*ones(1,hor);
           for ll = 1:hor
               tmp1 = (tmp_v(1,ll)^2+tmp_v(2,ll)^2)-desDist(i,t)^2;
%                tmp2 = tmp_vpeak(1,ll)^2+tmp_vpeak(2,ll)^2;
               if ll == hor
                   tmp_obj = tmp_obj+tmp1^2;%+0.3*tmp2^2;
               else
                   tmp_obj = tmp_obj+tmp1^2;%+0.3*tmp2^2;
               end
           end
%            tmp_obj = sum((sum(tmp_v.^2,1)-desDist(i,t)^2).^2);
%            tmp_obj = tmp_obj + 0*sum(sum((x(:,2:end)-peak*ones(1,hor)).^2));
%            for ll = 1:size(tmp_v,2)
%                tmp_obj = tmp_obj + tmp_v(1,ll)^2+tmp_v(2,ll)^2-2*sqrtm(tmp_v(1,ll)^2+tmp_v(2,ll)^2)*d(i,t)+d(i,t)^2;
%            end
        end
        obj3 = tmp_obj;
        %}
        
        %% obj4: regulation term to minimize the control input
        obj4 = 0;
        for j = 1:hor
            obj4 = obj4+u(1,j)^2+u(2,j)^2;
        end
        
        %% Define obj and constraints
        % define obj
        
        obj_w = [0.3;0.1;1;1];
        obj = obj_w'*[obj1;obj2;obj3;obj4];
        
        % define constraints
        constr = [x(:,1) == [rbt(i).x;rbt(i).y]];
        %
        for j = 1:hor
            % kinematic model
%             constr = [constr,x(1:2,j+1) == x(1:2,j)+u(1,j)*[cos(u(2,j));sin(u(2,j))]];
%             constr = [constr,vl <= u(1,j) <= vu];
            constr = [constr,x(:,j+1) == x(:,j)+u(:,j)];
            constr = [constr,x(1:2,j+1) >= [1;1] , x(1:2,j+1) <= [fld.x;fld.y]];
            constr = [constr,u(1,j)^2+u(2,j)^2 <= rbt(i).speedLimit^2];
        end
        %}
        
        % solve mpc
%         optset = sdpsettings('solver','fmincon','usex0',1,'debug',0,'verbose',0,...
%         'fmincon.Algorithm','sqp','fmincon.Display','final','fmincon.Diagnostics','off',...
%         'fmincon.TolCon',1e-5,'fmincon.TolFun',1e-5,'fmincon.MaxFunEvals',5000);
        optset = sdpsettings('solver','snopt','usex0',0,'debug',1,'verbose',1);
        sol = optimize(constr,obj,optset);
        if sol.problem == 0
            opt_x = value(x);
            opt_u = value(u);
        else
            err.time = [err.time,count];
            err.rbt = [err.rbt,i];
%             error(sprintf('fail to solve mpc for robot %d',i))
%             opt_x = init_x;
%             opt_u = init_u;
        end                
        
        % predict one-step robot motion and send it to 
%         pre_x = opt_x(1:2,end)+u(1,end)*[cos(u(2,end));sin(u(2,end))];
        pre_x = opt_x(1:2,end)+opt_u(:,end);
        rbt(i).plan_path = [opt_x(:,2:end),pre_x];
        rbt(i).x = opt_x(1,2);
        rbt(i).y = opt_x(2,2);
        rbt(i).opt_x{count} = opt_x;
        rbt(i).opt_u{count} = opt_u;
        
%         % calculate moving direction in unity
%         direction.x=rowOfPeak - rbt(i).x;
%         direction.y=colOfPeak - rbt(i).y;
%         lengthXY=sqrt(direction.x^2 + direction.y^2);
%         direction.x=direction.x/lengthXY;
%         direction.y=direction.y/lengthXY;
%         % move to peak PDF
%         %         rbt(i).x = rbt(i).x + rbt(i).speed * direction.x;
%         %         rbt(i).y = rbt(i).y + rbt(i).speed * direction.y;
%         DistCheck=0;
%         for t=rbt(i).neighbour
%             DistWithNeigh = sqrt((rbt(i).x - rbt(t).x)^2 + (rbt(i).y - rbt(t).y)^2);
%             if DistWithNeigh < 5 % Distance with neighbours are too short
%                 DistCheck = DistCheck + 1;
%             end
%         end
%         if (DistCheck  <= 0) % Enough distance
%             rbt(i).x = rbt(i).x + rbt(i).speed * direction.x;
%             rbt(i).y = rbt(i).y + rbt(i).speed * direction.y;
%         else % Too short distance--> go to random direction
%             rbt(i).x = rbt(i).x - rbt(i).speed * (rand(1)-0.5);
%             rbt(i).y = rbt(i).y - rbt(i).speed * (rand(1)-0.5);
%         end
    end
    
    %% Plot Planned Path
    % Plot plan_map
    hf3 = figure(3); clf(hf3)
    for i=1:NumOfRobot     
        subplot(2,3,i); contourf((rbt(i).plan_map)'); hold on;
        title(['Sensor ',num2str(i), ' Observ.= ',num2str(rbt(i).z)]);
        for j=1:NumOfRobot
            if i==j
                plot(rbt(j).opt_x{count}(1,1), rbt(j).opt_x{count}(2,1), 's','Color',rbt(j).color,'MarkerSize',8,'LineWidth',3);
                plot(rbt(j).opt_x{count}(1,2:end), rbt(j).opt_x{count}(2,2:end), 's','Color',rbt(j).color,'MarkerSize',5,'LineWidth',3);
            else
                plot(rbt(j).opt_x{count}(1,1), rbt(j).opt_x{count}(2,1), '.','Color',rbt(j).color,'MarkerSize',4,'LineWidth',1.5);
                plot(rbt(j).opt_x{count}(1,2:end), rbt(j).opt_x{count}(2,2:end), 'p','Color',rbt(j).color,'MarkerSize',5,'LineWidth',1.5);
            end
            plot(fld.tx, fld.ty, 'c+','MarkerSize',8,'LineWidth',3);
        end
    end
    subplot(2,3,5); xlabel(['Step=',num2str(count)]);
%     pause
    
    % save plots
    %
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