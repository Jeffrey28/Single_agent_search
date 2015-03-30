function outPara = pathPlanner(agent,inPara)
% include IPOPT in YALMIP
% addpath('D:\Program Files\MATLAB\2013a_crack\IPOPT3.11.8');
% addpath('D:\Chang Liu\ipopt');
% define input arguments
% x_h = inPara.pre_traj; % predicted human trajectory
hor = inPara.hor;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
% h_v = inPara.h_v;
obs_info = inPara.obs_info;
safe_marg = inPara.safe_marg;
campus = inPara.campus;
x_h = inPara.pre_traj;
all_comb = inPara.all_comb;
k = inPara.k;
% max_pts = inPara.max_pts;
prob_map_pf = inPara.prob_map_pf;
guess_u = inPara.guess_u;
guess_x = inPara.guess_x;
n_gmm = inPara.n_gmm;
last_r_obj = inPara.last_r_obj;
last_obj_w = inPara.last_obj_w;
chg_clt_flag = inPara.chg_clt_flag;

% define parameters
non_intersect_flag = 0; % flag for showing whether imposing the non-intersection constraint
dt = 0.05; % time interval for sampling the points on the line of the robot's path. used for imposing non-intersection constriant
safe_marg2 = 0.1; % margin for the robot's path line from the obstacle
tmp_hor = hor;
tc_scale = 1e-4; % scale for terminal cost
% h_v_value = norm(h_v,2);

init_state = [agent.sigma_s\agent.currentPos(1:2);agent.currentPos(3)];
inPara_cg = struct('hor',tmp_hor,'mpc_dt',mpc_dt,'safe_dis',safe_dis,...
    'safe_marg',safe_marg,'obs_info',obs_info,...
    'non_intersect_flag',non_intersect_flag,...
    'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state,...
    'campus',campus,'tc_scale',tc_scale,'x_h',x_h,'all_comb',{all_comb},'k',k,...
    'prob_map_pf',prob_map_pf,'guess_u',guess_u,'guess_x',guess_x,...
    'n_gmm',n_gmm,'last_r_obj',last_r_obj,'last_obj_w',last_obj_w,...
    'chg_clt_flag',chg_clt_flag); %'x',x,'u',u,'max_pts',max_pts,
% generate obj and constraints. contain a parameter that decides whether
% using the non-intersection constraints
tic;
[new_state,opt_u,opt_obj,obj_w] = solveMPC(inPara_cg);
display('MPC takes time:')
toc

outPara = struct('new_state',new_state,'opt_u',opt_u,'opt_obj',opt_obj,'obj_w',obj_w);
end

function [new_state,opt_u,opt_obj,obj_w] = solveMPC(inPara)
hor = inPara.hor;
% x = inPara.x;
% u = inPara.u;
% h_v = inPara.h_v;
mpc_dt = inPara.mpc_dt;
safe_dis = inPara.safe_dis;
safe_marg = inPara.safe_marg;
x_h = inPara.x_h;
obs_info = inPara.obs_info;
non_intersect_flag = inPara.non_intersect_flag;
% obj = inPara.obj;
% constr = inPara.constr;
agent = inPara.agent;
dt = inPara.dt;
safe_marg2 = inPara.safe_marg2;
% intgr_step = inPara.intgr_step;
campus = inPara.campus;
% tc_scale = inPara.tc_scale;
init_state = inPara.init_state;
all_comb = inPara.all_comb;
k = inPara.k;
% max_pts = inPara.max_pts;
prob_map_pf = inPara.prob_map_pf;
guess_u = inPara.guess_u;
guess_x = inPara.guess_x;
n_gmm = inPara.n_gmm;
last_r_obj = inPara.last_r_obj;
last_obj_w = inPara.last_obj_w;
chg_clt_flag = inPara.chg_clt_flag;

cur_clt = agent.cur_clt;
% hp_pt = agent.hp_pt;
grid_step = campus.grid_step;
% init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);

%% objective function
opt_obj = -100; % initialize the optimal obj value
% compact format
w = campus.w;
lambda = campus.lambda;
psi = campus.psi;
% prob_map = agent.updateProbMap(campus);
k_s = agent.k_s;

% remove terms with very small weights to speed up the optimization
%{
max_w = max(w);
rv_id = (abs(w) < max_w/10);
w(rv_id) = [];
sprintf('number of w is %d',length(w));
w = w/sum(w);
lambda(:,rv_id) = [];
psi(:,:,rv_id) = [];
%}
persistent x u Af

if k == 1    
    x = sdpvar(3,hor+1); %[lambda(2-D),theta]
    u = sdpvar(2,hor); %[v,a]
    
    Af = sdpvar(hor,1); % A funcitons [A_1,...,A_hor], A_1,...A_hor are for future positions
    for ii = 1:hor
        Af(ii) = A_fct2s(agent,x(1:2,ii+1),agent.psi_s);
    end
end

% note: obj1 is the probability of NOT detecting targets, we want to
% minimize this term
% w1 = 1;
% obj1 = 0;

% when the robot is in low probability area, don't use pond as obj. use
% obj2,3,4 to move to high probability area.
if (last_r_obj(2) < 0.98)
    w1 = 1;
else
    w1 = 0;
end
obj1 =0;
for jj = 1:length(w)
    alpha = sdpvar(size(all_comb,1),1);
    Af1 = A_fct2s(agent,lambda(:,jj),psi(:,:,jj));
    tmp_obj = 0;
    tmp_psi_prod = psi(:,:,jj);
    for kk = 1:size(all_comb,1) % give the term in the form of exp(A_(ij))-exp(A_i)-exp(A_j)+1
        tmp_idx = all_comb{kk};
        tmp_x = x(1:2,tmp_idx+1); % get the corresponding x
        tmp_para = updPara2([lambda(:,jj),tmp_x],...
            cat(3,psi(:,:,jj),repmat(agent.psi_s,[1,1,length(tmp_idx)])));
        tmp_psi = tmp_para.psi;
        tmp_lambda = tmp_para.lambda;
%         for ll = 1:length(tmp_idx)
            tmp_psi_prod = tmp_psi_prod*(agent.psi_s)^length(tmp_idx);
%         end
        tmp_psi_prod = tmp_psi_prod/tmp_psi;
        % it seems that if the robot is stuck in a small region for a few
        % steps, the tmp_psi_prod tends to be singular. Haven't looked into
        % why this happens        
        if abs(det(tmp_psi_prod)) < 1e-6
            tmp_psi_prod = tmp_psi_prod+eye(2)*1e-6;
        end
        tmp_coef = sqrt(det(tmp_psi_prod));
        alpha(kk) = A_fct2s(agent,tmp_lambda,tmp_psi)-Af1-sum(Af(tmp_idx));
        tmp_obj = tmp_obj+(-1)^length(tmp_idx)*k_s^length(tmp_idx)*tmp_coef*exp(alpha(kk));
    end
    if is(tmp_obj,'complex')
       display('complex obj') 
    end
    obj1 = obj1+w(jj)*(1+tmp_obj);
end
    
% add artificial potential function to the obstacles (static ones and
% humans)
w2 = 1;
obj2 = 0;
% potential function
%{
pen = 100; % rescale the gaussian distribution to impose penalty
pf_sigma_h = ((safe_dis+agent.d))^2*eye(2); % covariance matrix for the potential field for the human
for ii = 1:hor
    obj2 = obj2 + pen*mvnpdf((agent.sigma_s*x(1:2,ii+1))',x_h(2:3,ii+1)',pf_sigma_h);
    for jj = 1:size(obs_info,2)
        pf_sigma_obs = ((obs_info(3,jj)+safe_marg+agent.d))^2*eye(2); % covariance matrix for the potential field for the obstacles
        obj2 = obj2+pen*mvnpdf((agent.sigma_s*x(1:2,ii+1))',obs_info(1:2,jj)',pf_sigma_obs);
    end
end
%} 

% barrier function
for ii = 1:hor
    tmp_vec1 = agent.sigma_s*x(1:2,ii+1)-x_h(2:3,ii+1);
    obj2 = obj2 -log (sum(tmp_vec1.*tmp_vec1)-(safe_dis+agent.d)^2);
    for jj = 1:size(obs_info,2)
        tmp_vec2 = agent.sigma_s*x(1:2,ii+1)-obs_info(1:2,jj);
        obj2 = obj2 -log (sum(tmp_vec2.*tmp_vec2)-(obs_info(3,jj)+safe_marg+agent.d)^2);
    end
end
obj2 = 0.01*obj2; % rescale this term so that that negative term does not affect much when the robot is far from other objects

% determine whether the robot is already in the cur_clt region.
% clt_res = prob_map_pf(4,:);
% tmp_pt = prob_map_pf(1:2,clt_res == cur_clt);

%{
clt_res = agent.clt_res(3,:);
tmp_pt = agent.clt_res(1:2,clt_res == cur_clt);
tmp_vec = tmp_pt - agent.currentPos(1:2)*ones(1,size(tmp_pt,2));
tmp_dis = sqrt(sum(tmp_vec.*tmp_vec,1));
tmp_min = min(tmp_dis);
tmp_min_pt = tmp_pt(:,tmp_dis == tmp_min); % tmp_min_pt may contain several points
% tmp_min_pt = [45;35];
sprintf('the nearest point in current cluster is [%d,%d]',tmp_min_pt(1,1),tmp_min_pt(2,1))
% if the agent is not in the cur_clt region, add a terminal cost in the obj.
w3 = last_obj_w(3);
if tmp_min < agent.currentV * mpc_dt
    w3 = 0;
    obj3 = 0; 
else
    tmp_vec2 = agent.sigma_s*x(1:2,end) - (tmp_min_pt(:,1));%-1; *grid_step
    obj3 = norm(tmp_vec2);
    if  w3 == 0
        w3 = 1;
    else
        w3 = w3*1.2;
    end    
end
%}
% for debug use
% w3 = 1;
% tmp_vec2 = agent.sigma_s*x(1:2,end) - (tmp_min_pt(:,1));
% obj3 = norm(tmp_vec2);
obj3 = 0;
w3 = 0;

% add a terminal cost that guides the robot to the nearest highest probability
% point in current cluster or the goal cluster
% tmp_vec3 = max_pts*grid_step - ones(size(max_pts,1),1)*agent.currentPos(1:2)';
% tmp_dis2 = sqrt(sum(tmp_vec3.*tmp_vec3,2));
% tmp_min2 = min(tmp_dis2);
% tmp_min_pt2 = max_pts(tmp_dis2 == tmp_min2,:);

% this snippet may not work well since the particle number cannot
% represent the real probability, may need density estimation
clt_res = prob_map_pf(4,:);
cur_prob_map_pf = prob_map_pf(:,clt_res == agent.cur_clt);
cur_max_cnt = max(cur_prob_map_pf(3,:));
cur_max_idx = cur_prob_map_pf(3,:) == cur_max_cnt;
sprintf('max count of particles is %d at k = %d',cur_max_cnt,k)
sprintf('the position is [%d,%d]',cur_prob_map_pf(1,cur_max_idx),cur_prob_map_pf(2,cur_max_idx))
%         tmp_idx = (prob_map_pf(1:2,agent.clt_res == agent.cur_clt)==max(tmp_cnt));
max_pts = cur_prob_map_pf(1:2,cur_max_idx);
tmp_dif = max_pts-agent.currentPos(1:2)*ones(1,size(max_pts,2));
% find the closest point
tmp_idx = sum(tmp_dif.^tmp_dif,1) == min(sum(tmp_dif.^tmp_dif,1));
tmp_min_pt2 = max_pts(:,tmp_idx);
% tmp_min_pt2 = [30;30];
tmp_vec4 = agent.sigma_s*x(1:2,end) - tmp_min_pt2(:,1);%-1)*grid_step

% add this terminal cost only when the POND is greater than a threshold
% adjust the weight: if in previous step the POND is greater than 0.98, then increase the w4
w4 = last_obj_w(4);

% for debug use
% w4 = 1;
% obj4 = norm(tmp_vec4);

if (last_r_obj(2) < 0.98) && (chg_clt_flag == 0)
    w4 = 0;
    obj4 =0;
else
    if w4 == 0
        w4 = 1;
    else
        w4 = w4*1.1;
    end
    obj4 = norm(tmp_vec4);%*1e-1
end

obj_w = [w1;w2;w3;w4];
n_obj_w = obj_w/sum(obj_w); % normalized obj_w
obj = n_obj_w'*[obj1;obj2;obj3;obj4];

%% constraints
% get initial guess for NLP
% guess_u = [agent.currentV;0]*ones(1,hor);
% guess_state = agent.updState2(agent.currentPos,guess_u,mpc_dt);
% guess_x = [agent.sigma_s\guess_state(1:2,:);guess_state(3,:)];
% guess_x = []; % initial solution for the nonlinear problem
% guess_u = [];

while (1)
    % impose constraints
    % initial condition
    constr = [x(:,1) == init_state];
    for ii = 1:hor
        % constraints on robot dynamics
        % nonlinear constraint
        constr = [constr,agent.sigma_s*x(1:2,ii+1) == agent.sigma_s*x(1:2,ii)+u(1,ii)*[cos(x(3,ii));sin(x(3,ii))]*mpc_dt,...
            x(3,ii+1) == x(3,ii) + u(2,ii)*mpc_dt,...
            agent.minV<=u(1,ii)<=agent.maxV,agent.w_lb<=u(2,ii)<=agent.w_ub];
        
        % linear constraint
%         constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c,...
%             agent.minV<=u(1,ii)<=agent.maxV,agent.w_lb<=u(2,ii)<=agent.w_ub,...
%             x(1:2,ii+1) >= 0];
    end
    
    % solve MPC
    if ~isempty(guess_x)
        assign(x,guess_x); % assign initial solution
        assign(u,guess_u); % assign initial input
    end
    opt = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',0);
    opt.Algorithm = 'interior-point';
    sol = optimize(constr,obj,opt);
    
    if sol.problem ~= 0
        % if the orignial problem does not have a solution, warn it fails
        display('Fail to solve the original MPC')
        sol.info
        yalmiperror(sol.problem)
        figure % use empty figure to show that this if statement is executed
        
        % follow the MATLAB's suggestion: find a feasible point and use as
        % the initial solution for the original problem
        tmp_obj1 = n_obj_w(2)*obj2;
        assign(x,guess_x); % assign initial solution
        assign(u,guess_u); % assign initial input
        tmp_opt1 = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',0);
        %         tmp_opt.Algorithm = 'interior-point';
        sol = optimize(constr,tmp_obj1,tmp_opt1);
        if sol.problem ~= 0
            % if cannot find the feasible point using the nonlinear model, 
            % use the linear model.
            display('Fail to solve the original MPC using NLP, will use NLP with 0 obj')
            
            tmp_obj2 = 0;
            assign(x,guess_x); % assign initial solution
            assign(u,guess_u); % assign initial input
            tmp_opt2 = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',0);
            sol = optimize(constr,tmp_obj2,tmp_opt2);
            
            if sol.problem ~= 0
                display('Fail to solve the original MPC using NLP, will use LP with 0 obj')
                % linearize system
                [A,B,c] = agent.linearize_model3([agent.currentV;agent.currentPos(3)],mpc_dt);
                tmp_constr3 = [x(:,1) == init_state];
                tmp_constr3 = [tmp_constr3,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c,...
                    agent.minV<=u(1,ii)<=agent.maxV,agent.w_lb<=u(2,ii)<=agent.w_ub,...
                    x(1:2,ii+1) >= 0];
                tmp_obj3 = 0;
                assign(x,guess_x); % assign initial solution
                assign(u,guess_u); % assign initial input
                tmp_opt3 = sdpsettings('solver','linprog','debug',1,'verbose',0);
                sol = optimize(tmp_constr3,tmp_obj3,tmp_opt3);
                
                if sol.problem ~= 0
                    % break and the program will
                    % error in agentMove.m due to the returned empty new_state
                    new_state = [];
                    opt_u = [];
                    break
                elseif sol.problem == 0
                    guess_x = value(x);
                    guess_u = value(u);
                end
            elseif sol.problem == 0
                guess_x = value(x);
                guess_u = value(u);
            end
        elseif sol.problem == 0
            guess_x = value(x);
            guess_u = value(u);
        end
    elseif sol.problem == 0
        % test whether the solution without collision avoidance constraints
        % is not colliding
        %         opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        opt_u
        opt_obj = [value(obj);value(obj1);value(obj2);value(obj3);value(obj4)];
        %         break
        % update robot state
        new_state = agent.updState2(agent.currentPos,...
            opt_u,mpc_dt); % contains current and future states
        break
    end
end
end

function [outPara] = updPara2(lambda,psi)
% This function is similar to updateProbPara. However, this one will deal
% with the case for several parameters (such as updating mu using mu_i, 
% mu_k, mu_k+1, etc...)
outPara.psi = sum(psi,3);
outPara.lambda = sum(lambda,2);
end

function cst_chk = collisionCheck(agent,x,x_h,obs_info,safe_dis,safe_marg,safe_marg2)
len = size(x,2);
cst_chk = 1;
for ii = 1:len-1
    if sum((x(1:2,ii+1)-x_h(1:2,ii+1)).^2) < (safe_dis+agent.d)^2
        cst_chk = 0;
        return
    end
    for jj = 1:size(obs_info,2)
        if sum((x(1:2,ii+1)-obs_info(1:2,jj)).^2) < (obs_info(3,jj)+safe_marg+agent.d)^2
            cst_chk = 0;
            return
        end
    end
end
end

function [outPara] = updPara(mu,sigma)
% This function is similar to updateProbPara. However, this one will deal
% with the case for several parameters (such as updating mu using mu_i, 
% mu_k, mu_k+1, etc...)
num = size(mu,2);
tmp1 = 0;
tmp2 = 0;
for ii = 1:num
    tmp = eye(2)/sigma(:,:,ii);
    tmp1 = tmp1+tmp*mu(:,ii);
    tmp2 = tmp2+tmp;
end
outPara.sigma = eye(2)/tmp2;
outPara.mu = outPara.sigma*tmp1;
end
function cost = termCost(inPara)
% use the distance to the nearest high probability positions as the
% terminal cost. need to rescale this cost.
prob_map = inPara.prob_map;
campus = inPara.campus;
x_r = inPara.x_r;
tc_scale = inPara.tc_scale;

xMin = campus.endpoints(1);
xMax = campus.endpoints(2);
yMin = campus.endpoints(3);
yMax = campus.endpoints(4);
step = campus.grid_step;

% find the positions with top n largest probabilities
n = 5;
tmp_value = sort(prob_map(:),'descend');
[x_idx,y_idx] = find(prob_map >= tmp_value(n));
x_axis = xMin+step/2:step:xMax-step/2;
y_axis = yMin+step/2:step:yMax-step/2;
dif = x_r*ones(1,length(x_idx)) - [x_axis(x_idx);y_axis(y_idx)];
dif = reshape(dif,numel(dif),1);
cost = min(dif'*dif)*tc_scale;
end

function [obj,constr] = genMPCinfea(inPara)
% solve optimization problem when the initial condition is infeasible
hor = inPara.hor;
x = inPara.x;
u = inPara.u;
h_v = inPara.h_v;
mpc_dt = inPara.mpc_dt;
safe_dis = inPara.safe_dis;
safe_marg = inPara.safe_marg;
x_h = inPara.x_h;
obs_info = inPara.obs_info;
non_intersect_flag = inPara.non_intersect_flag;
obj = inPara.obj;
constr = inPara.constr;
agent = inPara.agent;
dt = inPara.dt;
safe_marg2 = inPara.safe_marg2;
init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);
% constr = [constr,sum((x(1:2,t)-x_h(:,t)).^2) >= safe_dis^2];
obj = -sum((x(1:2,2)-x_h(:,2)).^2);
for ii = 1:hor
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii) + u(1,ii)*mpc_dt, x(4,ii+1) == x(4,ii)+u(2,ii)*mpc_dt,...
       x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,agent.w_lb<=u(1,ii)<=agent.w_ub];
%     constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c];
%     constr = [constr,x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,-agent.maxW<=u(1,ii)<=agent.maxW];%
    % constraint on safe distance
    
%     constr = [constr,max(x(1:2,ii+1)-x_h(:,ii+1)) >= safe_dis];
    % constraint on obstacle avoidance
    % robot should not be inside the obstacles, i.e. robot waypoints should
    % not be inside the obstacle and the line connecting the waypoints 
    % should not intersect with the obstacle
%     [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
%{
    for jj = 1:size(obs_info,2)
        % waypoints not inside the obstacle
        constr = [constr,sum((x(1:2,ii+1)-obs_info(1:2,jj)).^2) >= (obs_info(3,jj)+safe_marg)^2];
        if non_intersect_flag == 1
            % line not intersecting with the obstacle
            n = floor(mpc_dt/dt);
            x0 = obs_info(1,jj); y0 = obs_info(2,jj);
            r = obs_info(3,jj);
            for kk = 0:n
                constr = [constr,sum((kk/n*x(1:2,ii+1)+(n-kk)/n*x(1:2,ii)-obs_info(1:2,jj)).^2)>=(r+safe_marg2)^2];
            end
        end
    end    
%}
end
end

function dis = dis_pl(p,v1,v2) % point-to-line distance
% define line equation: ax+by+c = 0
a = v2(2)-v1(2);
b = v1(1)-v2(1);
c = v1(2)*v2(1)-v2(2)*v1(1);
% poin-to-line distance
dis = (a*p(1)+b*p(2)+c)^2/(a^2+b^2);
end

function [a,b,c] = getLine(v1,v2)
% define line equation: ax+by+c = 0
a = v2(2)-v1(2);
b = v1(1)-v2(1);
c = v1(2)*v2(1)-v2(2)*v1(1);
end