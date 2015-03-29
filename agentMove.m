% 11/24/14
% simulate the movement of each agent
function [outPara] = agentMove(inPara)
%% initialization
% get input arguments
campus = inPara.campus;
agents = inPara.agents;
% way_pts = inPara.way_pts;

obv_traj = inPara.obv_traj;
est_state = inPara.est_state;
pre_traj = inPara.pre_traj;
plan_state = inPara.plan_state;
r_state = inPara.r_state;
r_input = inPara.r_input;
k = inPara.k;
hor = inPara.hor;
pre_type = inPara.pre_type;
samp_rate = inPara.samp_rate;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
safe_marg = inPara.safe_marg;
agentIndex = inPara.agentIndex;
plan_type = inPara.plan_type;



%% agents move 
% for agentIndex = 1:length(agents)
    agent = agents(agentIndex);
    
    %% human moves
    if strcmp(agent.type,'human') 
        h_tar_wp = inPara.h_tar_wp; % get human current target waypoint
        
        h_next_actions = getNextActionWithFixedHeading(agent.currentPos...
            ,h_tar_wp,agent.currentV,0,mpc_dt); % last argument is for zero deg_dev
        
        if k == 1
            tmp_agent_traj = agent.currentPos;
        else
            tmp_agent_traj = agent.traj;
        end
        
        % generate observations
        cur_pos = agent.currentPos(1:2);
        %
        % compute human heading
        %{
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
             tmp_t = (k-1)*mpc_dt+ii/samp_rate;
             obv_traj = [obv_traj,[tmp_t;cur_pos+(next_pos-cur_pos)*ii/samp_num]];
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
        %
        % estimated current state
%         est_state(:,k) = obv_traj(2:3,(k-1)*samp_rate*mpc_dt+1); 
%             est_state([1,3],k) = obv_traj(2:3,(k-1)*samp_rate*mpc_dt+1);
%             hd = cur_hd;
%             est_state([2,4],k) = h.currentV*[cos(hd);sin(hd)];
        %}
        
        %%  predict human future path
        %
        % prediction by GP
        samp_num = inPara.samp_num; % not used in current code. just put here so the outPara will not cause an error.
        pre_cov = inPara.pre_cov;
        r_obj = inPara.r_obj;
        if strcmp(pre_type,'GP')
            inPara_phj = struct('obv_traj',obv_traj,'hor',hor,'pre_type',pre_type,...
                'mpc_dt',mpc_dt,'samp_rate',samp_rate,'pre_cov',pre_cov);
            outPara_phj = predictHumanTraj(agent,inPara_phj);
            pre_traj(:,:,k) = outPara_phj.pre_traj;
            pre_cov(:,:,:,k) = outPara_phj.pre_cov;
%             pre_traj(:,:,k) = [[x_est((k-1)*samp_num+1,1);y_est((k-1)*samp_num+1,1)],[x_pre(k,:);y_pre(k,:)]];
%             pre_traj(:,:,k) = [x_pos_pre_imm(:,k)';y_pos_pre_imm(:,k)'];

        % prediction by extrapolation
        elseif strcmp(pre_type,'extpol')
%             inPara_phj = struct('state',est_state(:,k),'hor',hor,'pre_type',pre_type,...
%                 'mpc_dt',mpc_dt);
            inPara_phj = struct('obv_traj',obv_traj,...
                'hor',hor,'pre_type',pre_type,'mpc_dt',mpc_dt);
            outPara_phj = predictHumanTraj(agent,inPara_phj);
            pre_traj(:,:,k) = outPara_phj.pre_traj;
        end     
%         pos_pre_imm = inPara.pos_pre_imm;
        %}
        %% robot path planning
        all_comb = inPara.all_comb;
        clt_num = inPara.clt_num;
        guess_u = inPara.guess_u;
        guess_x = inPara.guess_x;
        n_gmm = inPara.n_gmm;
        obj_w = inPara.obj_w;
        % record current trajectory before moving
%         if k == 1
%             tmp_agent_traj = agent.currentPos;
%         else
            tmp_agent_traj = agent.traj;
%         end
        % clustering data
        %{
        prob_map = inPara.prob_map;
        clt_thresh = inPara.clt_thresh;
        if k == 1
            [agent.clt_res,agent.hp_pt] = agent.mapCluster(prob_map,clt_thresh);
        end
        %}
        prob_map_pf = inPara.prob_map_pf;
        if k == 1
            % [agent.clt_res,agent.hp_pt] = agent.mapCluster(prob_map_pf);
            [agent.clt_res,prob_map_pf] = agent.mapCluster(prob_map_pf,clt_num);
        end
        
        % decide which cluster to go. current strategy: go to the nearest
        % one
        %{
        agent.cur_clt = selectCluster(agent,campus.grid_step,prob_map);
        %}
        agent.cur_clt = selectCluster(agent,campus,prob_map_pf);
        
        % find the maximum probability point in current cluster
        %{
        cur_clt_pt = agent.hp_pt(agent.clt_res==agent.cur_clt,:);
        cur_clt_pt_idx = sub2ind(size(prob_map),cur_clt_pt(:,1),cur_clt_pt(:,2));
        cur_clt_prob = prob_map(cur_clt_pt_idx);
        [max_x,max_y] = ind2sub(size(prob_map),cur_clt_pt_idx(cur_clt_prob == max(cur_clt_prob)));
        %}
        % this snippet may not work well since the particle number cannot
        % represent the real probability, may need density estimation
%         clt_res = prob_map_pf(4,:);
%         cur_prob_map_pf = prob_map_pf(:,clt_res == agent.cur_clt);
%         [cur_max_cnt,cur_max_idx] = max(cur_prob_map_pf(3,:));
%         sprintf('max count of particles is %d at k = %d',cur_max_cnt(1),k)
%         %         tmp_idx = (prob_map_pf(1:2,agent.clt_res == agent.cur_clt)==max(tmp_cnt));
%         max_pts = cur_prob_map_pf(1:2,cur_max_idx);
        
        if strcmp(plan_type,'mpc')
%             inPara_pp = struct('pre_traj',pos_pre_imm(:,:,k),'hor',hor,...
%                 'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',[x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)],...
%                 'obs_info',campus.obs_info,'safe_marg',safe_marg);
            agent.currentPos = r_state(:,k); % update robot position and orientation
            
            if k > 1
                agent.currentV = r_input(1,k-1); % update robot speed
            end
            if k > 1
                last_r_obj = r_obj(:,k-1);
                last_obj_w = obj_w(:,k-1);
            else
                last_r_obj = zeros(size(r_obj,1),1);
                last_obj_w = zeros(size(obj_w,1),1);
            end
            
            inPara_pp = struct('hor',hor,'mpc_dt',mpc_dt,'campus',campus,...
                'obs_info',campus.obs_info,'safe_marg',safe_marg,'pre_traj',pre_traj(:,:,k),...
                'safe_dis',safe_dis,'all_comb',{all_comb},'k',k,...
                'prob_map_pf',prob_map_pf,'guess_u',guess_u,'guess_x',guess_x,...
                'n_gmm',n_gmm,'last_r_obj',last_r_obj,'last_obj_w',last_obj_w);%,[max_x,max_y],'max_pts',max_pts
            outPara_pp = pathPlanner(agent,inPara_pp);
%             opt_x = outPara_pp.opt_x;
            new_state = outPara_pp.new_state;
            opt_u = outPara_pp.opt_u;
            opt_obj = outPara_pp.opt_obj;
            obj_w(:,k) = outPara_pp.obj_w;
%             new_state = agent.updState([agent.currentPos(1:2);agent.currentV;agent.currentPos(3)],...
%                 opt_u,mpc_dt); % contains current and future states
%             agent.currentPos = [new_state(1:2,2);new_state(4,2)]; % update robot position and orientation
%             agent.currentV = new_state(3,2); % update robot speed
            r_state(:,k+1) = new_state(:,2); % save the robot's next state
            r_input(:,k) = opt_u(:,1);
            plan_state(:,:,k) = new_state;
            r_obj(:,k) = opt_obj;
            
            
            guess_u = [opt_u(:,2:end),zeros(size(opt_u,1),1)];
            guess_x = [new_state(:,2:end),new_state(:,end)];
            idx1 = (prob_map_pf(4,:) == 1);
            idx2 = (prob_map_pf(4,:) == 2);
            sprintf('# of particles in cluster 1 is %d',sum(prob_map_pf(3,idx1)))
            sprintf('# of particles in cluster 2 is %d',sum(prob_map_pf(3,idx2)))
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
            r_state(:,k+1) = opt_x(:,2);
            r_input(:,k) = opt_u(:,1);
            plan_state(:,:,k) = opt_x;
            
        end
        %}
        
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
            
    else
        error('Invalid agent type for planning path')
    end
% end
%% define output arguments
outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
    'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
    'samp_num',samp_num);
if exist('pre_cov', 'var')
    outPara.pre_cov = pre_cov;
end
if exist('r_obj', 'var')
    outPara.r_obj = r_obj;
end    
if exist('guess_u', 'var')
    outPara.guess_u = guess_u;
end  
if exist('guess_x', 'var')
    outPara.guess_x = guess_x;
end  
if exist('obj_w', 'var')
    outPara.obj_w = obj_w;
end  
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