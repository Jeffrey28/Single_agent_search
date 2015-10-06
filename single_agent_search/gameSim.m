% 11/24/2014
% this file is for the ME 290J course project

% 3/2/2015
% adjust the file for human-robot search scenario

clc 
clear % clear global variables
close all

%% Setup
% addpath('D:\Chang Liu\ipopt');
addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles');
addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
scale = 0.5; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');% docked

%%% define agents %%%
% setting for one-cluster
%{
% Human agent 1
h = agent('human');
h.currentPos = [22;20;0];%*scale% [x y heading]
% h.maxV = 1.5;
h.currentV = 2;

% Robot agent 1
r = agent('robot');
r.currentPos = [25;15;pi/2];%*scale; %[23.5;0.5;0];
r.currentV = 1;
r.maxV = 3;
r.minV = 0;
% r.a_lb = -1; 
% r.a_ub = 1;
r.w_lb = -pi;
r.w_ub = pi;
r.sigma_s = 9*eye(2);
r.psi_s = 1/2*eye(2)/r.sigma_s;
r.k_s = 2*pi*sqrt(det(r.sigma_s));
r.cur_clt = 0; % current goal cluster
r.d = 0.5; %robot diameter 

% define field type
campus_type = 'one_clt'; % one cluster campus
tar_pos = [35;40];
%}

% setting for two-cluster
%
% Human agent 1
h = agent('human');
h.currentPos = [22;20;0];%*scale% [x y heading]
% h.maxV = 1.5;
h.currentV = 2;

% Robot agent 1
r = agent('robot');
r.currentPos = [25;15;pi/2];%*scale; %[23.5;0.5;0];
r.currentV = 1;
r.maxV = 3;
r.minV = 0;
% r.a_lb = -1; 
% r.a_ub = 1;
r.w_lb = -pi/2;
r.w_ub = pi/2;
r.sigma_s = 9*eye(2);
r.psi_s = 1/2*eye(2)/r.sigma_s;
r.k_s = 2*pi*sqrt(det(r.sigma_s));
r.cur_clt = 0; % current goal cluster
r.d = 0.5; %robot diameter 

% define field type
campus_type = 'one_clt'; % one cluster campus
tar_pos = [27;40];
%}

%%% Set field %%%
xLength = 50;%*scale; 
yLength = 50;%*scale; 
xMin = 0;
yMin = 0;
xMax = xMin+xLength;
yMax = yMin+yLength;
grid_step = 0.5; % the side length of a probability cell
n_data = 3000;% number of randomly generated particles
inPara_gc = struct('campus_type',campus_type,'scale',scale,'n_data',n_data,...
    'xMin',xMin,'xMax',xMax,'yMin',yMin,'yMax',yMax,'grid_step',grid_step,...
    'tar_pos',tar_pos,'xLength',xLength,'yLength',yLength);
[campus,particles] = genCampus(inPara_gc);

% draw agents on the initial environment
figure(1)
hold on 
h1 = plot(h.currentPos(1),h.currentPos(2),'r','markers',2);
set(h1,'MarkerFaceColor','r');
set(h1,'MarkerEdgeColor','r');
set(h1,'Color','r');
set(h1,'LineStyle','-');
set(h1,'Marker','o');
h2 = plot(r.currentPos(1),r.currentPos(2),'g','markers',2);
set(h2,'MarkerFaceColor','g');
set(h2,'MarkerEdgeColor','g');
set(h2,'Color','g');
set(h2,'LineStyle','-');
set(h2,'Marker','o');

% generate obstacle points for drawing plots
c_set = campus.obs_info(1:2,:);
r_set = campus.obs_info(3,:);
theta_set = {{0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi}};
inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',theta_set,'type','obs');
obs_pos = getWayPts(inPara_gwp);

% draw obstacles on the initial environment
for jj = 1:size(obs_pos)
    tmp_pos = obs_pos{jj};
    fill(tmp_pos(1,:),tmp_pos(2,:),'b');
end
%}
xlim([0,campus.endpoints(2)]);
ylim([0,campus.endpoints(4)]);

%% Simulation
% simulation parameters
kf = 100; % simulation length (/s)
agents = [h,r];
hor = 3; % MPC horizon 
pre_type = 'extpol';%'extpol'; % 'extpol','IMM'. specify the method for predicting human motion
plan_type = 'mpc'; % 'MPC','greedy1','greedy0'. specify the method for robot controller
samp_rate = 20; % sampling rate (/Hz)
safe_dis = 0.5; %safe distance between human and robot
safe_marg = 0; % safety margin between human the the obstacle
mpc_dt = 0.5; % sampling time for model discretization used in MPC
prob_thresh = 0.8;
n_gmm = 7; % max cluster number

% define cluster number
switch campus_type
    case 'one_clt'
        clt_num = 1; % clustering number
        h.currentPos = [22;20;0];
        r.currentPos = [25;15;pi/2];
    case 'two_clt'
        clt_num = 2;
        h.currentPos = [20;30;0];
        r.currentPos = [18;7;pi/2];
end

% precompute combinatorial matrix
all_comb = {};
for ii = 1:hor
    all_comb = [all_comb;num2cell(nchoosek(1:hor,ii),2)];
end

% initialize variables
% obv_traj = zeros(3,0); % observed human trajectory; first row denotes the time. [t,x,y]
plan_traj = zeros(3,hor+1,kf); % predicted trajectory of neighboring agents [t,x,y]
plan_state = zeros(3,hor+1,kf); % robot's current and planned future state [x,y,v]
r_state = zeros(3,kf); % robot's actual state [x,y,theta]
r_state(:,1) = r.currentPos;
r_input = zeros(2,kf); % robot's actual control input [v,a]
r_obj = zeros(5,kf); % save the objective funciton values [obj;obj1;obj2;obj3;obj4]
wp_cnt = 1; % the waypoint that the human is heading for
h_tar_wp = h_way_pts(:,wp_cnt); % the way point that the human is heading for
sensor_reading = -1*ones(length(agents),kf); % record the sensor readings of each agent
prob_map_set = []; % record probability map for each stepR
tar_found = 0; % binary variable. 1 indicates that the target is found
clt_thresh = 1e-4; %threshold for choosing points to be clustered
pre_cov = zeros(hor,hor,hor,kf); % covariance of the predicted x and y trajectory.
F = []; % save frame of plots for movie
% particles = zeros(2,n_data); % save the generated particles;
guess_x = []; % guess of the initial solution of the MPC
guess_u = []; % guess of the initial solution of the MPC
obj_w = zeros(4,kf); % save the weights for the objective functions
% addpath('.\sim_res')
% load('x_pos_pre_imm','x_pos_pre_imm')
% load('y_pos_pre_imm','y_pos_pre_imm')
% pos_pre_imm = zeros(2,size(x_pos_pre_imm,1)-1,size(x_pos_pre_imm,2));
% for ii = 1:size(x_pos_pre_imm,2)
%     pos_pre_imm(:,:,ii) = [x_pos_pre_imm(2:end,ii)';y_pos_pre_imm(2:end,ii)'];
% end
for k = 1:kf
    display(k)
    %% generate particles
%     if k == 1
%        particles = agents(2).genParticles(campus.w,campus.mu,campus.sigma,n_data); 
%     end
    %% update probability map
    % Set prior probabilities based on agent positions and target readings
    w = campus.w;
    mu = campus.mu;
    sigma = campus.sigma;
    
    % Generate sensor readings
    for agentIndex = 1:length(agents)
        agent = agents(agentIndex);
        if (strcmp(agent.type,'robot'))
            cur_pos = agent.currentPos(1:2);
            [~,~,sim_reading] = sensorSim(agent,tar_pos,cur_pos);
            sensor_reading(agentIndex,k) = sim_reading;
            sprintf('sensor reading is %d',sim_reading)
            % update parameters for probability map using update law
            %{
            cur_lambda = agent.sigma_s\cur_pos;
            cur_psi = 1/2*eye(2)/agent.sigma_s;
            [w,mu,sigma] = agent.updateProbPara(w,lambda,psi,sim_reading,cur_pos,agent.sigma_s);
            inPara_upp = struct('w',w,'lambda',lambda,'psi',psi,'obs',sim_reading,...
                'lambda_s',cur_lambda,'psi_s',cur_psi,'field',campus);
            outPara_upp = agent.updateProbPara2(inPara_upp);
            w = outPara_upp.w;
            mu = outPara_upp.mu;
            sigma = outPara_upp.sigma;
            lambda = outPara_upp.lambda;
            psi = outPara_upp.psi; 
            old_w = outPara_upp.old_w;
            old_mu = outPara_upp.old_mu;
            old_sigma = outPara_upp.old_sigma;
            old_lambda = outPara_upp.old_lambda;
            old_psi = outPara_upp.old_psi; 
            %}
            % update probability map using particle filter
            inPara_pf = struct('particles',particles,'obs',sim_reading,'n_gmm',n_gmm);
            outPara_pf = agent.particle_filter(inPara_pf);
            particles = outPara_pf.particles;
            % following quantities are from the fitted gmm
            campus.w = outPara_pf.w;
            campus.lambda = outPara_pf.lambda;
            campus.psi = outPara_pf.psi;
            campus.sigma = outPara_pf.sigma;
            campus.mu = outPara_pf.mu;
        end
    end
    
    % draw the pf-based probability map
    figure
    plot(particles(1,:),particles(2,:),'o')
    
    % remove terms with very small weights
    %{
    [max_w,max_id] = max(w);
    rv_id = (abs(w) < max(w)*1e-3);
    w(rv_id) = [];
    w = w/sum(w);
    lambda(:,rv_id) = [];
    psi(:,:,rv_id) = [];
    mu(:,rv_id) = [];
    sigma(:,:,rv_id) = [];
    %}
    % this part is used with the updateProbPara2 to save the original map
    % parameters based on update law and the new ones based on gmm fitting
    %{
    campus.w = w;
    campus.lambda = lambda;
    campus.psi = psi;
    campus.sigma = sigma;
    campus.mu = mu;
    tmp_campus = campus;
    tmp_campus.w = old_w;
    tmp_campus.lambda = old_lambda;
    tmp_campus.psi = old_psi;
    tmp_campus.sigma = old_sigma;
    tmp_campus.mu = old_mu;
    %}
    
    % get the two probability distribution: the fitted one and the one
    % before fitting
%     [tmp_prob_map,~] = agent.updateProbMap(tmp_campus,0);
    [prob_map_pf,prob_map_gmm] = agent.updateProbMap_pf(campus,particles);
  %     prob_map_set(:,:,k) = matrixToCartesian(prob_map);
    
    %% Test to see if game is over
    
    % Game over if an agent is on the same grid as the target and recieves a
    % detection. Human agent will be the first agent in the array agents.
    %{
    for ii = 1:length(agents)
        agent = agents(ii);
%         if strcmp(agent.type,'robot') && (norm((agent.currentPos(1:2)-campus.targetPos),1) <= agent.currentV*mpc_dt ...
%             && sim_reading(ii) == 1)
%             sprintf('Target has been found! Game ends at t=%d.',t)
%             tar_found = 1;
%             break       
%         end
        for jj = 1:length(max_id)
            if norm(sigma(:,:,max_id(jj))) <= 3
                sprintf('Target has been found! Game ends at k=%d.',k)
                tar_found = 1;
                return       
            end
        end
    end
    %}
    
    inPara_ec = struct('prob_thresh',prob_thresh,'prob_map_pf',prob_map_pf,...
        'r',agents(2),'campus',campus,'sensor_reading',sensor_reading(2,1:k),...
        'particles',particles);
    game_end = endCheck(inPara_ec);
    if game_end == 1
        sprintf('Target has been found! Game ends at k=%d.',k)
        return
    end
    
    %% search
    if tar_found == 0
        %% change human speed
        %{
    tmp_v = agents(1).currentV+h_acl(k)*mpc_dt;
    if tmp_v <= 0
        tmp_v = 0;
    else
        agents(1).currentV = tmp_v;
    end
        %}
        
        %% waypoint check
        % check if the human needs to choose the next way point
        %
        inPara_ec = struct('h',agents(1),'h_way_pts',h_way_pts,'wp_cnt',wp_cnt,'mpc_dt',mpc_dt);
        [outPara_ec] = wpCheck(inPara_ec);
        
        wp_end = outPara_ec.wp_end;
        arr_wp_flag = outPara_ec.arr_wp_flag; % the human has arrived at a way point
        
        if wp_end == 1
            sprintf('human has arrived at his destinaiton at k = %d',k-1) 
        else            
            if (k == 1)
                h_tar_wp = h_way_pts(:,1);               
            end
            if arr_wp_flag == 1
                wp_cnt = wp_cnt+1;
                h_tar_wp = h_way_pts(:,wp_cnt);
                %         agents(1).currentV = h_v(wp_cnt);
            end
        end
        %}
        %% both agents move
        % human moves and robot predicts and plans its own path
        
        % human moves
        %
        agentIndex = 1;
        inPara_ams = struct('campus',campus,'agents',agents,...
            'obv_traj',obv_traj,'est_state',est_state,...
            'pre_traj',plan_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
            'k',k,'hor',hor,'pre_type',pre_type,'samp_rate',samp_rate,...
            'safe_dis',safe_dis,'mpc_dt',mpc_dt,'safe_marg',safe_marg,...
            'agentIndex',agentIndex,'plan_type',plan_type,'h_tar_wp',h_tar_wp);
        [outPara_ams] = agentMove(inPara_ams);
        agents = outPara_ams.agents;
        obv_traj = outPara_ams.obv_traj;
        samp_num = outPara_ams.samp_num;
        %}
        
        % robot moves
        %
        agentIndex = 2;
        %     load('obv_traj3_w_time.mat')% Load Tracjectory of Human
        %     obv_traj1=obv_traj';
        %     parameter;  % parameter for IMM
        inPara_ams = struct('campus',campus,'agents',agents,'h_tar_wp',h_tar_wp,...
            'obv_traj',obv_traj,'est_state',est_state,...
            'pre_traj',plan_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
            'k',k,'hor',hor,'pre_type',pre_type,'samp_rate',samp_rate,...
            'safe_dis',safe_dis,'mpc_dt',mpc_dt,'safe_marg',safe_marg,...
            'agentIndex',agentIndex,'plan_type',plan_type,'samp_num',samp_num,...
            'prob_map_pf',prob_map_pf,'clt_thresh',clt_thresh,'pre_cov',pre_cov,...
            'all_comb',{all_comb},'clt_num',clt_num,'r_obj',r_obj,...
            'guess_u',guess_u,'guess_x',guess_x,'n_gmm',n_gmm,'obj_w',obj_w);
        
        [outPara_ams] = agentMove(inPara_ams);
        agents = outPara_ams.agents;
        est_state = outPara_ams.est_state;
        plan_traj = outPara_ams.pre_traj;
        pre_cov = outPara_ams.pre_cov;
        plan_state = outPara_ams.plan_state;
        r_state = outPara_ams.r_state;
        r_input = outPara_ams.r_input;
        r_obj = outPara_ams.r_obj;
        guess_x = outPara_ams.guess_x;
        guess_u = outPara_ams.guess_u;
        obj_w = outPara_ams.obj_w;
        disp('plan_state is:')
        disp(plan_state(:,:,k))
        disp('r_obj is:')
        disp(r_obj(:,k))
        disp('obj_w is:')
        disp(obj_w(:,k))
        %}
    end
    
    %% plot trajectories
    
    % plot specifications
    color_agent = {'r','g','r','g'};
    marker_agent = {'o','o','^','^'};
    line_agent = {'-','-',':',':'};
    line_w_agent = [3,3,3,3];
    orange = [1 204/255 0];
    color_target = {'m','b',orange};
    figure;
%     drawnow
    hold on
    box
    grid on
    axis equal
    
    % draw probability map
    %{
    plot_prob_map = [prob_map zeros(size(prob_map,1),1); zeros(1,size(prob_map,2)) 0]';
    p_handle = pcolor(xMin:grid_step:xMax,yMin:grid_step:yMax,plot_prob_map);
    set(p_handle,'EdgeColor','none');
    colormap(b2r_m(min(plot_prob_map(:)),max(plot_prob_map(:))));
    colorbar
    %}
    
    % draw gmm probability map
    plot_prob_map_gmm = [prob_map_gmm zeros(size(prob_map_gmm,1),1); zeros(1,size(prob_map_gmm,2)) 0]';
    p_handle = pcolor(xMin:grid_step:xMax,yMin:grid_step:yMax,plot_prob_map_gmm);
    set(p_handle,'EdgeColor','none');
    colormap(b2r_m(min(plot_prob_map_gmm(:)),max(plot_prob_map_gmm(:))));
%     colormap(b2r_m(0,1));
%     caxis([min(plot_prob_map_gmm(:)),max(plot_prob_map_gmm(:))]);
%     colormap('default');
    colorbar
    
    % draw targets
    %
    for jj = 1:campus.targetNum
        h = plot(campus.targetPos(1,jj),campus.targetPos(2,jj),'MarkerSize',30);
        set(h,'Marker','p');
    end
    %}
    
    % draw obstacles
    %
    for jj = 1:size(obs_pos)
        tmp_pos = obs_pos{jj};
        fill(tmp_pos(1,:),tmp_pos(2,:),'b');
    end
    %}
    % make the plot square
    xlim([0,campus.endpoints(2)]);
    ylim([0,campus.endpoints(4)]);
    
    % draw agent trajectory
    %
    
    for ii = 1:length(agents)
        tmp_agent = agents(ii);
        h1 = plot(tmp_agent.traj(1,:),tmp_agent.traj(2,:),'markers',2);
        set(h1,'MarkerFaceColor',color_agent{ii});
        set(h1,'MarkerEdgeColor',color_agent{ii});
        set(h1,'Color',color_agent{ii});
        set(h1,'LineStyle',line_agent{ii});
        set(h1,'Marker',marker_agent{ii});
        set(h1,'LineWidth',line_w_agent(ii));
        h2 = plot(tmp_agent.currentPos(1),tmp_agent.currentPos(2),color_agent{ii},'markers',2);
        set(h2,'MarkerFaceColor',color_agent{ii});
        set(h2,'MarkerEdgeColor',color_agent{ii});
        set(h2,'Color',color_agent{ii});
        set(h2,'LineStyle',line_agent{ii});
        set(h2,'Marker',marker_agent{ii});
        % get hte points to draw the safe region around the human and the size
        % of the robot
        % human
%         if strcmp(tmp_agent.type,'human')
%             c_set = [tmp_agent.traj(1,:);tmp_agent.traj(2,:)];
%             r_set = safe_dis;
%             theta = 0:pi/8:2*pi;
%             inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta',theta,'type','agent');
%             circle_pos = getWayPts(inPara_gwp);
%             for jj = 1:size(circle_pos)
%                 tmp_pos = circle_pos{jj};
%                 fill(tmp_pos(1,:),tmp_pos(2,:),color_agent{ii});
%             end
%         elseif strcmp(tmp_agent.type,'robot')
            % the following line can be problematic, need to change after
            % fixing the problem in time series of sensor data
%             c_set = [tmp_agent.traj(1,:);tmp_agent.traj(2,:)];
%             r_set = tmp_agent.d;
%             theta = 0:pi/8:2*pi;
%             inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta',theta,'type','agent');
%             circle_pos = getWayPts(inPara_gwp);
%             for jj = 1:size(circle_pos)
%                 tmp_pos = circle_pos{jj};
%                 fill(tmp_pos(1,:),tmp_pos(2,:),color_agent{ii});
%             end
%         end
    end
    %}
    
    % draw the sensor range (sig)
    r = agents(2);
    c_set = r.currentPos(1:2);
    r_set = sqrt(r.sigma_s(1,1));
    theta = 0:pi/8:2*pi;
    inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta',theta,'type','agent');
    circle_pos = getWayPts(inPara_gwp);
    for jj = 1:size(circle_pos)
        tmp_pos = circle_pos{jj};
        plot(tmp_pos(1,:),tmp_pos(2,:),'--b');
    end
    
    % predicted human positions
    %
    h3 = plot(plan_traj(2,:,k),plan_traj(3,:,k),color_agent{3},'markers',2);
    set(h3,'MarkerFaceColor',color_agent{3});
    set(h3,'MarkerEdgeColor',color_agent{3});
    set(h3,'Color',color_agent{3});
    set(h3,'LineStyle',line_agent{3});
%     set(h3,'Marker',marker_agent{3});
    set(h3,'LineWidth',3);%line_w_agent(3)
    c_set = [plan_traj(2,2:end,k);plan_traj(3,2:end,k)];
    r_set = safe_dis;
    theta = 0:pi/8:2*pi;
    inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta',theta,'type','agent');
    circle_pos = getWayPts(inPara_gwp);
    for jj = 1:size(circle_pos)
        tmp_pos = circle_pos{jj};
        plot(tmp_pos(1,:),tmp_pos(2,:));%,color_agent{3});
    end
    %}
    
    % planned robot trajectory
    %
    if (tar_found == 0)
        
        h4 = plot(plan_state(1,:,k),plan_state(2,:,k),color_agent{2},'markers',2);
        set(h4,'MarkerFaceColor',color_agent{4});
        set(h4,'MarkerEdgeColor',color_agent{4});
        set(h4,'Color',color_agent{4});
        set(h4,'LineStyle',line_agent{4});
%         set(h4,'Marker',marker_agent{4});
        set(h4,'LineWidth',3);%line_w_agent(4)
        c_set = [plan_state(1,2:end,k);plan_state(2,2:end,k)];
        r_set = r.d;
        theta = 0:pi/8:2*pi;
        inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta',theta,'type','agent');
        circle_pos = getWayPts(inPara_gwp);
        for jj = 1:size(circle_pos)
            tmp_pos = circle_pos{jj};
            plot(tmp_pos(1,:),tmp_pos(2,:));%,color_agent{4});
        end
        grid minor
        set(gca,'MinorGridLineStyle','-','XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5])
        axis equal
        xlim([0,xLength]);ylim([0,yLength]);
    end
    % close plots when there are too many plots
    h5 = gcf;
    if h5 > 30
        close all;
    end
    %}
    
    % draw the difference map between the fitted and original maps
    %{
    figure
    dif_prob_map = prob_map-tmp_prob_map;
    if sum(dif_prob_map(:)) ~= 0
        
        plot_dif_prob_map = [dif_prob_map zeros(size(prob_map,1),1); zeros(1,size(prob_map,2)) 0]';
        p_handle = pcolor(xMin:grid_step:xMax,yMin:grid_step:yMax,plot_dif_prob_map);
        set(p_handle,'EdgeColor','none');
        colormap(b2r_m(min(plot_dif_prob_map(:)),max(plot_dif_prob_map(:))));
        colorbar
    end
    %}
    
    
    
    % draw animation of the search process. does not work well
    %{
    x_fig = 450;
    y_fig = 150;
    w_fig = 650;% width
    h_fig = 600;% height
    hd_fig = gcf;
    set(hd_fig, 'Position', [x_fig y_fig w_fig h_fig])
    F = [F,getframe(hd_fig)];
    %}
   
end

% pre_traj = pos_pre_imm;
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