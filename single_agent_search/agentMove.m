% 11/24/14
% simulate the movement of each agent
function [outPara] = agentMove(inPara)
%% initialization
% get input arguments
campus = inPara.campus;
agents = inPara.agents;
plan_traj = inPara.plan_traj; % planned trajectory of neighboring agents
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

%% agent move 

%% update probability map using map and planned path from neighboring agents
% combine the probability map

% update the probability map based on neighboring agents' planned path


agent = agents(agentIndex);

%%  predict human future path
%{
% prediction by GP
samp_num = inPara.samp_num; % not used in current code. just put here so the outPara will not cause an error.
pre_cov = inPara.pre_cov;
r_obj = inPara.r_obj;
if strcmp(pre_type,'GP')
    inPara_phj = struct('obv_traj',obv_traj,'hor',hor,'pre_type',pre_type,...
        'mpc_dt',mpc_dt,'samp_rate',samp_rate,'pre_cov',pre_cov);
    outPara_phj = predictHumanTraj(agent,inPara_phj);
    plan_traj(:,:,k) = outPara_phj.pre_traj;
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
    plan_traj(:,:,k) = outPara_phj.pre_traj;
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
tmp_agent_traj = agent.traj;

% clustering data
prob_map_pf = inPara.prob_map_pf;
if k == 1
    [agent.clt_res,prob_map_pf] = agent.mapCluster(prob_map_pf,clt_num);
end

% decide which cluster to go. current strategy: go to the nearest one
[agent.cur_clt,chg_clt_flag] = selectCluster(agent,campus,prob_map_pf);

if strcmp(plan_type,'mpc')
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
        'obs_info',campus.obs_info,'safe_marg',safe_marg,'pre_traj',plan_traj(:,:,k),...
        'safe_dis',safe_dis,'all_comb',{all_comb},'k',k,...
        'prob_map_pf',prob_map_pf,'guess_u',guess_u,'guess_x',guess_x,...
        'n_gmm',n_gmm,'last_r_obj',last_r_obj,'last_obj_w',last_obj_w,...
        'chg_clt_flag',chg_clt_flag);
    outPara_pp = pathPlssanner(agent,inPara_pp);
    new_state = outPara_pp.new_state;
    opt_u = outPara_pp.opt_u;
    opt_obj = outPara_pp.opt_obj;
    obj_w(:,k) = outPara_pp.obj_w;
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
% end
%% define output arguments
outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
    'pre_traj',plan_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
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