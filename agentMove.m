% 11/24/14
% simulate the movement of each agent
function [outPara] = agentMove(inPara)
%% initialization
% get input arguments
campus = inPara.campus;
agents = inPara.agents;
% way_pts = inPara.way_pts;
% h_tar_wp = inPara.h_tar_wp;
% obv_traj = inPara.obv_traj;
% est_state = inPara.est_state;
% pre_traj = inPara.pre_traj;
plan_state = inPara.plan_state;
r_state = inPara.r_state;
r_input = inPara.r_input;
k = inPara.k;
hor = inPara.hor;
% pre_type = inPara.pre_type;
% samp_rate = inPara.samp_rate;
% safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
safe_marg = inPara.safe_marg;
agentIndex = inPara.agentIndex;
plan_type = inPara.plan_type;

%% agents move 
% for agentIndex = 1:length(agents)
    agent = agents(agentIndex);
    
    %% human moves
    if strcmp(agent.type,'human')         
        h_next_actions = getNextActionWithFixedHeading(agent.currentPos...
            ,h_tar_wp,agent.currentV,0,mpc_dt); % last argument is for zero deg_dev
        
        if k == 1
            tmp_agent_traj = agent.currentPos;
%             obv_traj = [0;agent.currentPos(1:2)];
        else
            tmp_agent_traj = agent.traj;
        end
        
        % generate observations
        cur_pos = agent.currentPos(1:2);
        %
        % compute human heading
        if k == 1
            % we assume that human heading at the beginning is known.
            cur_hd = agent.currentPos(3)+h_next_actions(3);
        else
            cur_hd = agent.currentPos(3);
        end
        %}
        agent = takeNextAction(agent,h_next_actions);
        next_pos = agent.currentPos(1:2);
        t = norm(cur_pos-next_pos,2)/agent.currentV; % calculate the time for human to move to his next position
        samp_num = double(uint16(t*samp_rate)); % get the number of observations of human position
        for ii = 1:samp_num
             % observed human position
             tmp_t = (k-1)*mpc_dt+(ii-1)/samp_rate;
             obv_traj = [obv_traj,[tmp_t;cur_pos+(next_pos-cur_pos)*(ii-1)/samp_num]];
        end
        
        % update human position
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
    %% robot predicts and moves
    elseif strcmp(agent.type,'robot')  
%         samp_num = inPara.samp_num;
%         h = agents(1);
%         cur_hd = h.currentPos(3);
        %% estimate human position
        %{
        [x_est,y_est,x_pre,y_pre,x_pos_est,input,time] = IMM_Com_run();
        est_state([1,2],:,k) = x_est((k-1)*samp_num+1:end-1,:)';
        est_state([3,4],:,k) = y_est((k-1)*samp_num+1:end-1,:)';
        %}
        % estimation with no measurement noise
        %{
            est_state([1,3],k) = obv_traj(2:3,(k-1)*samp_rate*mpc_dt+1);
            hd = cur_hd;
            est_state([2,4],k) = h.currentV*[cos(hd);sin(hd)];
        %}
        
        %%  predict human future path
        %{
        % prediction by IMM
        if strcmp(pre_type,'IMM')
            pre_traj(:,:,k) = [[x_est((k-1)*samp_num+1,1);y_est((k-1)*samp_num+1,1)],[x_pre(k,:);y_pre(k,:)]];
%             pre_traj(:,:,k) = [x_pos_pre_imm(:,k)';y_pos_pre_imm(:,k)'];
        % prediction by extrapolation
        elseif strcmp(pre_type,'extpol')
%             inPara_phj = struct('state',est_state(:,k),'hor',hor,'pre_type',pre_type,...
%                 'mpc_dt',mpc_dt);
            inPara_phj = struct('state',[x_est((k-1)*samp_num+1,:)';y_est((k-1)*samp_num+1,:)'],...
                'hor',hor,'pre_type',pre_type,'mpc_dt',mpc_dt);
            pre_traj(:,:,k) = predictHumanTraj(agent,inPara_phj);
        end     
        pos_pre_imm = inPara.pos_pre_imm;
        %}
        %% robot path planning
        %
        if strcmp(plan_type,'mpc')
%             inPara_pp = struct('pre_traj',pos_pre_imm(:,:,k),'hor',hor,...
%                 'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',[x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)],...
%                 'obs_info',campus.obs_info,'safe_marg',safe_marg);
            inPara_pp = struct('hor',hor,'mpc_dt',mpc_dt,'campus',campus,...
                'obs_info',campus.obs_info,'safe_marg',safe_marg);
            outPara_pp = pathPlanner(agent,inPara_pp);
            opt_x = outPara_pp.opt_x;
            opt_u = outPara_pp.opt_u;
            agent.currentPos = opt_x(1:3,2); % update robot position
            agent.currentV = opt_x(4,2); % update robot speed
            r_state(:,k) = opt_x(:,2);
            r_input(:,k) = opt_u(:,1);
            plan_state(:,:,k) = opt_x;
        elseif strcmp(plan_type,'greedy1') || strcmp(plan_type,'greedy0')
            inPara_pp = struct('pre_traj',pos_pre_imm(:,:,k),'hor',hor,...
                'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',...
                [x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)],'obs_info',campus.obs_info,...
                'safe_marg',safe_marg,'plan_type',plan_type);
            outPara_pp = pathPlannerGreedy(agent,inPara_pp);
            opt_x = outPara_pp.opt_x;
            opt_u = outPara_pp.opt_u;
            agent.currentPos = opt_x(1:3,2); % robot moves
            agent.currentV = opt_x(4,2); % robot updates its speed
            r_state(:,k) = opt_x(:,2);
            r_input(:,k) = opt_u(:,1);
            plan_state(:,:,k) = opt_x;
        end
        %}
        if k == 1
            tmp_agent_traj = agent.currentPos;
        else
            tmp_agent_traj = agent.traj;
        end
        
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
            
    else
        error('Invalid agent type for planning path')
    end
% end
%% define output arguments
% outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
%     'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
%     'samp_num',samp_num);
outPara = struct('agents',agents,'plan_state',plan_state,'r_state',r_state,...
    'r_input',r_input);
end

function next_act = getNextActionWithFixedHeading(a_pos,t_pos,v,deg_dev,mpc_dt)
% calculate the actions for agent to move from a_pos to target position
% t_pos with velocity v
vec = t_pos - a_pos(1:2);
heading = calAngle(vec)+deg_dev;
% next_acts = zeros(3); %[dx,dy,d_heading]. only the first action changes the agent's heading
tmp_a_pos = a_pos;
dx = v*cos(heading)*mpc_dt;
dy = v*sin(heading)*mpc_dt;
d_hd = heading - tmp_a_pos(3);
next_act = [dx;dy;d_hd];
end    