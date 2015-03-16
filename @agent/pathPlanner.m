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

% define parameters
non_intersect_flag = 0; % flag for showing whether imposing the non-intersection constraint
dt = 0.05; % time interval for sampling the points on the line of the robot's path. used for imposing non-intersection constriant
safe_marg2 = 0.1; % margin for the robot's path line from the obstacle
tmp_hor = hor;
tc_scale = 1e-4; % scale for terminal cost
% h_v_value = norm(h_v,2);

init_state = [agent.sigma_s\agent.currentPos(1:2);agent.currentV;agent.currentPos(3)];
% while(tmp_hor > 0)
    % define MPC
    x = sdpvar(4,tmp_hor+1); %[lambda(2-D),theta,v]
    u = sdpvar(2,tmp_hor); %[w,a]
    
    % constraints on future states
%     inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'h_v',h_v_value,'mpc_dt',mpc_dt,...
%         'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
%         'non_intersect_flag',non_intersect_flag,'obj',1,'constr',constr,...
%         'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state,...
%         'intgr_step',intgr_step);
    inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'mpc_dt',mpc_dt,'safe_dis',safe_dis,...
        'safe_marg',safe_marg,'obs_info',obs_info,...
        'non_intersect_flag',non_intersect_flag,...
        'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state,...
        'campus',campus,'tc_scale',tc_scale,'x_h',x_h);
    % generate obj and constraints. contain a parameter that decides whether
    % using the non-intersection constraints
    [new_state,opt_u] = solveMPC(inPara_cg); 
    
%     % solve MPC
%     opt = sdpsettings('solver','fmincon','usex0',1,'debug',1);
%     opt.Algorithm = 'interior-point';
%     sol = optimize(constr,obj,opt);
%     
%     if sol.problem == 0
%         opt_x = value(x); % current and future states
%         opt_u = value(u); % future input
% %         break
%     else
%         display('Fail to solve MPC')
%         sol.info
%         yalmiperror(sol.problem)
%         return
%     end
    %{
    if sol.problem == 0
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        if tmp_hor < hor
            opt_x = [opt_x,opt_x(:,end)*ones(1,hor-tmp_hor)]; % current and future states
            opt_u = [opt_u,zeros(size(opt_u,1),hor-tmp_hor)]; % future input
        end
        for ii = 1:tmp_hor
            for jj = 1:size(obs_info,2)
                % check if the line intersects with some obstacle
                n = floor(mpc_dt/dt);
                % x0 = obs_info(1,jj); y0 = obs_info(2,jj);
                r = obs_info(3,jj);
                for kk = 0:n
                    tmp = sum((kk/n*opt_x(1:2,ii+1)+(n-kk)/n*opt_x(1:2,ii)-obs_info(1:2,jj)).^2) - (r+safe_marg2)^2;
                    if tmp < 0
                        non_intersect_flag = 1;
                        break
                    end
                end
                if tmp < 0
                    break
                end
            end
            if tmp < 0
                break
            end
        end
        if tmp >= 0
            break
        end
    else
        display('Fail to solve MPC')
        sol.info
        yalmiperror(sol.problem)
        % safe_dis = safe_dis/2;   
        tmp_hor = tmp_hor-1;
    end
    %}
% end

%{
if tmp_hor == 0 % if the MPC fails, just find the input at the next step to maximize the humna-robot distance
    %{
    % define MPC
    x = sdpvar(4,hor+1); %[x,y,theta,v]
    u = sdpvar(2,hor); %[w,a]
    tmp_hor = 1;
    % impose constraints
    % initial condition
    constr = [x(:,1) == init_state];
    % constraints on future states
    %     for tmp_hor = 2:hor+1
    inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'h_v',h_v_value,'mpc_dt',mpc_dt,...
        'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
        'non_intersect_flag',non_intersect_flag,'obj',0,'constr',constr,...
        'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state);
    [obj,constr] = genMPCinfea(inPara_cg); % generate constraints. contain a parameter that decides whether using the non-intersection constraints
    % solve MPC
    opt = sdpsettings('solver','ipopt','usex0',1,'debug',1);
    sol = optimize(constr,obj,opt);
    
    if sol.problem == 0
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        opt_x(:,tmp_hor+2:end) = zeros(size(opt_x,1),size(opt_x,2)-tmp_hor-1);
        opt_u(:,tmp_hor+2:end) = zeros(size(opt_u,1),size(opt_u,2)-tmp_hor-1);
        %{
        for ii = 1:hor
            for jj = 1:size(obs_info,2)
                % check if the line intersects with some obstacle
                n = floor(mpc_dt/dt);
                %                 x0 = obs_info(1,jj); y0 = obs_info(2,jj);
                r = obs_info(3,jj);
                for kk = 0:n
                    tmp = sum((kk/n*opt_x(1:2,ii+1)+(n-kk)/n*opt_x(1:2,ii)-obs_info(1:2,jj)).^2) - (r+safe_marg2)^2;
                    if tmp < 0
                        non_intersect_flag = 1;
                        break
                    end
                end
                if tmp < 0
                    break
                end
            end
            if tmp < 0
                break
            end
        end
        if tmp >= 0
            break
        end
        %}
    else
        display('Fail to solve MPC')
        sol.info
        yalmiperror(sol.problem)
        %         safe_dis = safe_dis/2;
        %         tmp_hor = tmp_hor-1;
    end
    %     end
    %}
    x_r = agent.currentPos(1:2);
    r_hd = agent.currentPos(3);
    r_v = agent.currentV;
    a_lb = agent.a_lb;
    a_ub = agent.a_ub;
    w_lb = agent.w_lb;
    w_ub = agent.w_ub;
    x_r_next = x_r+r_v*[cos(r_hd);sin(r_hd)]*mpc_dt;
    rh_dis_next = sqrt(sum((x_r_next - x_h(:,2)).^2));
    rh_dir = calAngle(x_h(:,2)-x_r_next); % direction from robot to human
    if (rh_dis_next >= safe_dis)
        % if robot will be outside of the collision region, then turn its
        % heading toward the human's next position
        min_hd = r_hd + w_lb*mpc_dt;
        max_hd = r_hd + w_ub*mpc_dt;
        if rh_dir<=max_hd && rh_dir>=min_hd
            r_hd_next = rh_dir;
        else
            hd_dif_min = min(abs(min_hd-rh_dir),abs(2*pi-abs(min_hd-rh_dir)));
            hd_dif_max = min(abs(max_hd-rh_dir),abs(2*pi-abs(max_hd-rh_dir)));
            if hd_dif_min < hd_dif_max
                r_hd_next = min_hd;
            else
                r_hd_next = max_hd;
            end
        end
        
        %     r_v_next = r_v + a_ub*mpc_dt;
    else
        % if robot will be inside collision region, then turn its
        % heading against the human's next position
        min_hd = r_hd + w_lb*mpc_dt;
        max_hd = r_hd + w_ub*mpc_dt;
        op_rh_dir = -rh_dir;
        op_rh_dir = op_rh_dir - floor(op_rh_dir/(2*pi))*2*pi; % opposite direction
        if op_rh_dir<=max_hd && op_rh_dir>=min_hd
            r_hd_next = op_rh_dir;
        else
            hd_dif_min = min(abs(min_hd-rh_dir),abs(2*pi-abs(min_hd-rh_dir)));
            hd_dif_max = min(abs(max_hd-rh_dir),abs(2*pi-abs(max_hd-rh_dir)));
            if hd_dif_min < hd_dif_max
                r_hd_next = max_hd;
            else
                r_hd_next = min_hd;
            end
            %     elseif op_rh_dir<min_hd
            %         r_hd_next = max_hd;
            %     elseif op_rh_dir>max_hd
            %         r_hd_next = min_hd;
        end
        
        %     r_v_next = r_v + a_lb*mpc_dt;
    end
    tmp = r_hd_next;
    tmp = tmp - 2*pi*floor(tmp/(2*pi));
    r_hd_next = tmp;


    % chang robot speed to match human's current estimated speed
    min_v = r_v + a_lb*mpc_dt;
    max_v = r_v + a_ub*mpc_dt;
    
    if (rh_dis_next >= 2*safe_dis)
        r_v_next = max_v;
    elseif (rh_dis_next >= safe_dis) && (rh_dis_next < 2*safe_dis)
        if norm(h_v,2) >= max_v
            r_v_next = max_v;
        elseif norm(h_v,2) <= min_v
            r_v_next = min_v;
        else
            r_v_next = norm(h_v,2);
        end
    else
        r_v_next = min_v;
    end
    r_v_next = max(r_v_next,0);
    
    opt_x = [[x_r;r_hd;r_v],[x_r_next*ones(1,hor);r_hd_next*ones(1,hor);r_v_next,zeros(1,hor-1)]];
    opt_u = [[(r_hd_next-r_hd)/mpc_dt;(r_v_next-r_v)/mpc_dt],zeros(2,hor-1)];
end
%}

outPara = struct('new_state',new_state,'opt_u',opt_u);
end

function [new_state,opt_u] = solveMPC(inPara)
hor = inPara.hor;
x = inPara.x;
u = inPara.u;
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

cur_clt = agent.cur_clt;
clt_res = agent.clt_res;
hp_pt = agent.hp_pt;
grid_step = campus.grid_step;
% init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);
% objective function
% integral format
%{
xMin = campus.endpoints(1);
xMax = campus.endpoints(2);
yMin = campus.endpoints(3);
yMax = campus.endpoints(4);
x_p = xMin+intgr_step/2:intgr_step:xMax-intgr_step/2; % integral points on x axis
y_p = yMin+intgr_step/2:intgr_step:yMax-intgr_step/2;
pond_mat = sdpvar(length(x_p),length(y_p)); % matrix of Pr(ND)
sigma_s_inv = eye(2)/agent.sigma_s;
for x1 = 1:length(x_p)
    for x2 = 1:length(y_p)
        x_t = [x_p(x1);y_p(x2)];
        for ii = 1:hor            
            x_r = x(1:2,ii+1);
            pond = 1-exp(-1/2*(x_t-x_r)'*sigma_s_inv*(x_t-x_r));
            if ii == 1
                pond_mat(x1,x2) = pond;
            else
                pond_mat(x1,x2) = pond_mat(x1,x2)*pond;
            end
        end
    end
end
obj = sum(pond_mat(:));
%}
% compact format
% this is for horizon of 2
w = campus.w;
lambda = campus.lambda;
psi = campus.psi;
% prob_map = agent.updateProbMap(campus);
k_s = agent.k_s;

% remove terms with very small weights to speed up the optimization
max_w = max(w);
rv_id = (abs(w) < max_w/20);
w(rv_id) = [];
sprintf('number of w is %d',length(w));
w = w/sum(w);
lambda(:,rv_id) = [];
psi(:,:,rv_id) = [];

A2 = A_fct2(agent,x(1:2,2),agent.psi_s);
A3 = A_fct2(agent,x(1:2,3),agent.psi_s);
obj1 = 0;
for jj = 1:length(w)
    A1 = A_fct2(agent,lambda(:,jj),psi(:,:,jj));
    alpha = sdpvar(3,1);
    tmp_para = updPara2([lambda(:,jj),x(1:2,2),x(1:2,3)],...
        cat(3,psi(:,:,jj),agent.psi_s,agent.psi_s));
    tmp_psi = tmp_para.psi;
    tmp_lambda = tmp_para.lambda;
    alpha(3) = A_fct2(agent,tmp_lambda,tmp_psi)-A1-A2-A3;
    tmp_para = updPara2([lambda(:,jj),x(1:2,3)],cat(3,psi(:,:,jj),agent.psi_s));
    tmp_psi = tmp_para.psi;
    tmp_lambda = tmp_para.lambda;
    alpha(2) = A_fct2(agent,tmp_lambda,tmp_psi)-A1-A3;
    tmp_para = updPara2([lambda(:,jj),x(1:2,2)],cat(3,psi(:,:,jj),agent.psi_s));
    tmp_psi = tmp_para.psi;
    tmp_lambda = tmp_para.lambda;
    alpha(1) = A_fct2(agent,tmp_lambda,tmp_psi)-A1-A2;
    obj1 = obj1+w(jj)*(1-k_s*exp(alpha(1))-k_s*exp(alpha(2))+k_s^2*exp(alpha(3)));
end

% add artificial potential function to the obstacles (static ones and
% humans)
obj2 = 0;
pf_sigma_h = ((safe_dis+agent.d)/3)^2*eye(2); % covariance matrix for the potential field for the human
for ii = 1:hor
   obj2 = obj2 + mvnpdf((agent.sigma_s*x(1:2,ii+1))',x_h(2:3,ii+1)',pf_sigma_h);
   for jj = 1:size(obs_info,2)
       pf_sigma_obs = ((obs_info(3,jj)+safe_marg+agent.d)/3)^2*eye(2); % covariance matrix for the potential field for the obstacles
       obj2 = obj2+mvnpdf((agent.sigma_s*x(1:2,ii+1))',obs_info(1:2,jj)',pf_sigma_obs);
   end
end

% determine whether the robot is already in the cur_clt region.
tmp_pt = hp_pt(clt_res == cur_clt,:);
tmp_vec = tmp_pt*grid_step - ones(size(tmp_pt,1),1)*agent.currentPos(1:2)';
tmp_dis = sqrt(sum(tmp_vec.*tmp_vec,2));
tmp_min = min(tmp_dis);
tmp_min_pt = tmp_pt(tmp_dis == tmp_min,:); % tmp_min_pt may contain several points
% if the agent is not in the cur_clt region, add a terminal cost in the obj.
obj3 =0;
if tmp_min > agent.currentV * mpc_dt 
    tmp_vec2 = agent.sigma_s*x(1:2,end) - tmp_min_pt(1,:)';
    obj3 = norm(tmp_vec2)*1e-2;
end

obj = 1/3*obj1+1/3*obj2+1/3*obj3;
% inPara_tc = struct('prob_map',prob_map,'x_r',x(1:2,end),'tc_scale',tc_scale,...
%     'campus',campus);
% obj = obj + termCost(inPara_tc);
cst_flag = 0;
init_x = []; % initial solution for the nonlinear problem
while (1)
    % constraints
    % linearize system
    [A,B,c] = agent.linearize_model2([agent.currentV;agent.currentPos(3)],mpc_dt);
    % impose constraints
    % initial condition
    constr = [x(:,1) == init_state];
    for ii = 1:hor
        % constraints on robot dynamics
        % nonlinear constraint
        %     constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(3,ii)*[cos(x(4,ii));sin(x(4,ii))]*mpc_dt,...
        %         x(3,ii+1) == x(3,ii) + u(1,ii)*mpc_dt, x(4,ii+1) == x(4,ii)+u(2,ii)*mpc_dt,...
        %         x(3,ii+1)>=0,agent.a_lb<=u(1,ii)<=agent.a_ub,agent.w_lb<=u(2,ii)<=agent.w_ub];
        
        % linear constraint
        constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c...
            0<=x(3,ii+1)<=agent.maxV,agent.a_lb<=u(1,ii)<=agent.a_ub,agent.w_lb<=u(2,ii)<=agent.w_ub];
        %{
        if cst_flag == 1
            % if need to include constraints
            % constraint on safe distance
            constr = [constr,sum((agent.sigma_s*x(1:2,ii+1)-x_h(2:3,ii+1)).^2) >= (safe_dis+agent.d)^2];
            %     constr = [constr,max(x(1:2,ii+1)-x_h(:,ii+1)) >= safe_dis];
            
            % constraint on obstacle avoidance
            % robot should not be inside the obstacles, i.e. robot waypoints should
            % not be inside the obstacle and the line connecting the waypoints
            % should not intersect with the obstacle
            %     [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
            %
            for jj = 1:size(obs_info,2)
                % waypoints not inside the obstacle
                constr = [constr,sum((agent.sigma_s*x(1:2,ii+1)-obs_info(1:2,jj)).^2) >= (obs_info(3,jj)+safe_marg+agent.d)^2];
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
        end
        %}
    end
    
    % solve MPC
    if ~isempty(init_x)
        assign(x,init_x); % assign initial solution
    end
    opt = sdpsettings('solver','fmincon','usex0',1,'debug',1);
    opt.Algorithm = 'interior-point';
    sol = optimize(constr,obj,opt);
    
    if sol.problem ~= 0 && cst_flag == 0
        % if the orignial problem does not have a solution, warn it fails
        display('Fail to solve the original MPC')
        sol.info
        yalmiperror(sol.problem)
        figure % use empty figure to show that this if statement is executed
        
        % follow the MATLAB's suggestion: find a feasible point and use as
        % the initial solution for the original problem
        tmp_obj = 0;
        opt = sdpsettings('solver','linprog','usex0',1,'debug',1);
%         opt.Algorithm = 'interior-point';
        sol = optimize(constr,tmp_obj,opt);
        if sol.problem ~= 0
            % if cannot find the feasible point, break and the program will
            % error in agentMove.m due to the returned empty new_state
            % canno
            new_state = [];
            opt_u = [];
            break
        elseif sol.problem == 0
            init_x = value(x);
        end
    elseif sol.problem == 0 && cst_flag == 0
        % test whether the solution without collision avoidance constraints
        % is not colliding
%         opt_x = value(x); % current and future states
        opt_u = value(u); % future input
%         break
        % update robot state
        new_state = agent.updState([agent.currentPos(1:2);agent.currentV;agent.currentPos(3)],...
            opt_u,mpc_dt); % contains current and future states
        break
        % check if the new state will violate the collision constraints
        %{
        cst_chk = collisionCheck(agent,new_state,x_h,obs_info,safe_dis,safe_marg,safe_marg2);
        if cst_chk == 0
            % if violating constraints, re-do the MPC, considering
            % collision aovidance
            cst_flag = 1;
        else
            break
        end
        %}
    elseif sol.problem ~= 0 && cst_flag == 1
        % if solution cannot be found after considering the collsion
        % constraints, then use the control policy from the unconstrained
        % case
        display('cannot satisfy the collision avoidance constraints, use the input from the unconstraint case')
        break
    elseif sol.problem == 0 && cst_flag == 1
        % get a solution for MPC with collision avoidance constraints
        opt_u = value(u); % future input
        new_state = agent.updState([agent.currentPos(1:2);agent.currentV;agent.currentPos(3)],...
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