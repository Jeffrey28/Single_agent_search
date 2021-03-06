 classdef Robot
    properties
        % motion specs
        traj;
        state; % [x;y;theta;v]
        a_lb;
        a_ub;
        w_lb;
        w_ub;
        v_lb;
        v_ub;
        
        % sensor specs
        sensor_type;
        theta0; % sensor range in angle
        r; % sensor range
        % linear sensor
        C; % C matrix of observation
        % nonlinear sensor
        h; % y=h(x)
        del_h; % gradient of h
        R; % covariance for sensor model
        
        % sensor modeling
%         gam; % function handle for gamma
%         gam_den; % function handle for the denominator of gamma
%         gam_aprx; % function handle for the linearized gamma
%         p_aprx; % function handle for the element of the linearized covriance matrix
        alp1; % parameters in gam
        alp2; % parameters in gam
        alp3; % parameters in gam
        thr; % the threshold to avoid very large exponential function value
        
        % trust region
        tr_inc; % increment factor of trust region
        tr_dec; % decrement factor of trust region
        % penalty factor
        mu_inc;
        
        % observation
        y; % observation measurement
        
        % target
        target;
        
        % filtering
        % xKF
        est_pos; % estimated target position
        P; % estimation covariance
        est_pos_hist;
        P_hist;
        % GMM
        gmm_num; % # of gmm components
        max_gmm_num; % max # of gmm components for gmm fitting purpose
        wt; % weigths of gmm components
        % PF
        particles; % array of particle positions [x;y]
        gmm_mu; % array of mean [x;y]
        gmm_sigma; % cell of covariance
        
        % path planning
        mpc_hor;
        dt;
        optu;
        
        % configuration of optimization
        cfg;
        
        % performance metrics
        ml_pos;
        ent_pos;
    end
    
    methods
        function this = Robot(inPara)
            this.state = inPara.state;
            this.traj = inPara.state;
            this.a_lb = inPara.a_lb;
            this.a_ub = inPara.a_ub;
            this.w_lb = inPara.w_lb;
            this.w_ub = inPara.w_ub;
            this.v_lb = inPara.v_lb;
            this.v_ub = inPara.v_ub;
            
            % sensor specs
            this.R = inPara.R;
            this.h = inPara.h;
            this.del_h = inPara.del_h;
            this.theta0 = inPara.theta0;
            this.r = inPara.range;
            
            % gamma modeling
%             this.gam = inPara.gam; % function handle for gamma
%             this.gam_den = inPara.gam_den; % function handle for the denominator of gamma
%             this.gam_aprx = inPara.gam_aprx; % function handle for the linearized gamma
%             this.p_aprx = inPara.p_aprx; % function handle for the element of the linearized covriance matrix
            this.alp1 = inPara.alp1; % parameters in gam
            this.alp2 = inPara.alp2;
            this.alp3 = inPara.alp3;
            this.thr = inPara.thr;
            
            % trust region
            this.tr_inc = inPara.tr_inc;
            this.tr_dec = inPara.tr_dec;
            % penalty factor
            this.mu_inc = inPara.mu_inc;
            
            % target
            this.target = inPara.target;
            
            % filtering
            this.sensor_type = inPara.sensor_type;
            % xKF
            this.est_pos = inPara.est_pos;
            this.P = inPara.P;
            this.est_pos_hist = [];
            this.P_hist = [];
            % gmm
%             this.gmm_num = inPara.gmm_num;
%             this.wt = inPara.wt;
            % pf
            this.max_gmm_num = inPara.max_gmm_num;
            this.particles = inPara.particles;
            
            this.mpc_hor = inPara.mpc_hor;
            this.dt = inPara.dt;
            this.optu = [];
            
            this.cfg = inPara.cfg;
        end
        
        %% sensor modeling
        % approximate straightline edge of sensor FOV based on current
        % robot state
        function [a,b] = FOV(this,st)
            theta = st(3);
            x0 = st(1);
            y0 = st(2);
            alp1 = theta - this.theta0;
            alp2 = theta + this.theta0;
            a = [sin(alp1),-cos(alp1);-sin(alp2),cos(alp2)]; % [a1;a2]
            b = [x0*sin(alp1)-y0*cos(alp1);-x0*sin(alp2)+y0*cos(alp2)];%[b1;b2];
        end
        
        % determine if the target is in sensor FOV
        function flag = inFOV(this,tar_pos)
            [a,b] = this.FOV(this.state);
            flag = (a(1,:)*tar_pos-b(1) <= 0) && (a(2,:)*tar_pos-b(2) <= 0)...
                && (norm(tar_pos-this.state(1:2)) <= this.r);            
        end
        
        %% measurement generation
        % generate a random measurement
        function y = sensorGen(this,fld)
            tar_pos = fld.target.pos;
            
            % draw FOV and target position. for debug purpose
            this.drawFOV(this.state,fld,'cur')
            hold on
            plot(fld.target.pos(1),fld.target.pos(2),'b','markers',5,'Marker','*');
                        
            if strcmp(this.sensor_type,'rb')
                % range-bearing sensor
                if this.inFOV(tar_pos)
                    y = this.h(tar_pos,this.state(1:2))+(mvnrnd([0;0],this.R))';
                    display(mvnpdf(y,this.h(tar_pos,this.state(1:2)),this.R))
                else
                    y = [-100;-100];
                end
            elseif strcmp(this.sensor_type,'ran')
                if this.inFOV(tar_pos)
                    y = norm(tar_pos-this.state(1:2))+normrnd(0,this.R);
                else
                    y = -100;
                end
            elseif strcmp(this.sensor_type,'lin')
                if this.inFOV(tar_pos)
                    y = this.h(tar_pos,this.state(1:2))+(mvnrnd([0;0],this.R))';                    
                else
                    y = [-100;-100];
                end
            end
        end
        
        %% filtering
        function this = KF(this,fld)
            % measurement
            y = this.y;
            
            % target
            tar = fld.target;
            A = tar.A;
            B = tar.B;
            
%             f = tar.f;
%             del_f = tar.del_f;
%             A = del_f;
            
            Q = tar.Q;

            % current estimation
            x = this.est_pos;
            P = this.P;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            
            % prediction
            x_pred = A*x+B;            
            P_pred = A*P*A'+Q;
            
            % update
            C = del_h(x,this.state(1:2));
            if sum(y-[-100;-100]) ~= 0
                % if an observation is obtained
                K = P_pred*C'/(C*P_pred*C'+R);
                x_next = x_pred+K*(y-h(x_pred,this.state(1:2)));%C*x_pred-(x_pred-this.state(1:2))
                P_next = P_pred-K*C*P_pred;
            else
                x_next = x_pred;
                P_next = P_pred;
            end
            
            this.est_pos = x_next;
            this.P = P_next;
            this.est_pos_hist = [this.est_pos_hist,x_next];
            this.P_hist = [this.P_hist,P_next];
            % this makes KF compatible with cvxPlanner_scp
            this.gmm_num = 1; 
            this.wt = 1;
        end
        
        function this = GSF(this,fld)
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            % measurement
            y = this.y;
            
            % used for updating gmm component weights
            alp = ones(this.gmm_num,1);
            
            for ii = 1:this.gmm_num
                % current estimation
                P = this.P{ii};
                x = this.est_pos(:,ii);
                A = del_f(x);
                % prediction
                x_pred = f(x); %%% replace this one with new nonlinear model
                P_pred = A*P*A'+Q;
                
                % update
                % sensor model linearization
                C = del_h(x_pred);
                
                if sum(y-[-100;-100]) ~= 0
                    % if an observation is obtained
                    K = P_pred*C'/(C*P_pred*C'+R);
                    x_next = x_pred+K*(y-h(x_pred));
                    P_next = P_pred-K*C*P_pred;
                    alp(ii) = mvnpdf(y,h(x_pred),C*P_pred*C'+R);
                else
                    x_next = x_pred;
                    P_next = P_pred;
                end
                
                this.est_pos(:,ii) = x_next;
                this.P{ii} = P_next;
            end
            
            % update gmm component weight
            wt = this.wt.*alp;
            this.wt = wt/sum(wt);
            tar_mean = this.est_pos*this.wt;
            tar_cov = zeros(2,2);
            for jj = 1:this.gmm_num
                tar_cov = tar_cov+this.wt(jj)*this.P{jj};
            end
            this.est_pos_hist = [this.est_pos_hist,tar_mean];
            this.P_hist = [this.P_hist,tar_cov];
        end
        
        function this = PF(this,fld)
            % particle filter
            
            % target
            tar = fld.target;
            f = tar.f;
            % sensor
            h = this.h;
            % measurement
            y = this.y;
            
            particles = this.particles;
            
            %% particle filtering
            np = size(particles,2); % number of particles
            
            % initalize particles weights
            w = zeros(np,1);
            
            % state update: since we use the static target, no update is needed
            pred_par = zeros(2,np); % predicted particle state
            for ii = 1:np
                pred_par(:,ii) = f(particles(:,ii));
            end
            
            % weight update
            for ii = 1:np
                if sum(y == -100) >= 1
                    % if the target is outside FOV.
                    if this.inFOV(pred_par(:,ii))
                        w(ii) = 10^-20;
                    else
                        w(ii) = 1;
                    end
                else
                    if this.inFOV(pred_par(:,ii))
                        w(ii) = mvnpdf(y,this.h(pred_par(:,ii),this.state(1:2)),this.R);
                    else
                        w(ii) = 10^-20;
                    end
                end
            end
            w = w/sum(w);
            
            % resampling
            idx = randsample(1:np, np, true, w);
            new_particles = particles(:,idx);
            this.particles = new_particles;
            
            %%%%% may need to deal with practical issues in PF, to be
            %%%%% filled later
            
            %% gmm fitting
            max_gmm_num = this.max_gmm_num; % maximum gmm component number
            
            gmm_model = cell(max_gmm_num,1);
            opt = statset('MaxIter',1000);
            AIC = zeros(max_gmm_num,1);
            tic;
            for kk = 1:max_gmm_num
                gmm_model{kk} = fitgmdist(new_particles',kk,'Options',opt,...
                    'Regularize',0.001,'CovType','full');
                AIC(kk)= gmm_model{kk}.AIC;
            end
            sprintf('Robot.m, Line %d. gmm fitting takes time as:',MFileLineNr());
            toc;
            
            [minAIC,numComponents] = min(AIC);
            
            best_model = gmm_model{numComponents};
            this.gmm_num = numComponents;
            this.gmm_mu = best_model.mu';
            this.gmm_sigma = best_model.Sigma;
            this.wt = best_model.PComponents';
            
            % convert the data form to be compatible with the main code
            this.est_pos = this.gmm_mu(:);
            for ii = 1:numComponents
                this.P{ii} = this.gmm_sigma(:,:,ii);
            end
        end
        
        %% planning
        
        % try re-writing the problem using cvx solver. Different from
        % cvxPlanner below, which formulates the problem as a convex P (turns
        % out not!), this one formulates the problem as QP each iteration
        function [optz,optu, merit, model_merit, new_merit] = cvxPlanner_scp(this,fld,optz,optu,plan_mode) % cvxPlanner(this,fld,init_sol)
            % merit, model_merit, new_merit are the merit function value of
            % x, approx merit value of xp, and merit function value of xp
            
            % use the multi-layer approach similar to Sachin's work. Fix
            % the parameter for the sensor, solve path planning. Then
            % refine the parameter until close to reality. In each
            % iteration, a QP program is solved. The initial solution
            % comes from ngPlanner
            
            % planing in non-Gaussian (GMM) belief space
            N = this.mpc_hor;
            dt = this.dt;
            
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            alp_inc = 5; % increament paramter for alphatr_inc
%             
            % config for optimization
            cfg = this.cfg;
            
            % set up simulation            
            if isempty(optz)
                prev_state = [];
            else
                prev_state = struct('optz',optz,'optu',optu);
            end
            
            switch plan_mode
                case 'lin'
                    init_sol = genInitState_kf(this,fld,prev_state);
                case 'nl'
                    init_sol = genInitState(this,fld,prev_state);
            end
            
            zref = init_sol.z;
            uref = init_sol.u;
            xref = init_sol.x;
            x_pred_ref = init_sol.x_pred;
            Kref = init_sol.K;
            Pref = init_sol.P;
            P_pred_ref = init_sol.P_pred;
            % construct the state
            % s is the collection of all decision variables
            % s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:)]   
            [s,snum] = this.setState(zref,uref,xref,x_pred_ref,Pref,P_pred_ref,Kref);
            
            cfg.snum = snum;
            
            %%% functions and parameters for scp                     
            %%% objective
%             obj = @(s) this.getObj(s,snum);
            
            % Q and q in quad linear obj term: x'*Q*x/2+q'*x
            fQ = zeros(length(s)); 
            q = zeros(1,length(s));            
            
            %%% kinematics constraints
            objKin = @(s) this.getKinConstr(s,snum); 
            
            %%% linear equality constraints
            %{
            A_eq = zeros(snum.total); 
            % initial condition on states
            % z(:,1) == this.state;
            % x(:,1) == this.est_pos(:);
            % for jj = 1:this.gmm_num
              %  triu(P(:,:,jj,1)) == triu(this.P{jj});
            % end            
            tmp_idx = [1:4,snum.idx(2,2)+1:snum.idx(2,2)+2*this.gmm_num,...
                snum.idx(4,2)+1:snum.idx(4,2)+4*this.gmm_num];
            for ii = tmp_idx
                A_eq(ii,ii) = 1;
            end
            vec = [];
            for jj = 1:this.gmm_num
                vec = [vec;this.P{jj}(:)];
            end
            b_eq = [this.state;zeros(snum.idx(2,2)-length(this.state),1);...
                this.est_pos(:);zeros(snum.idx(4,2)-snum.idx(2,2)-length(this.est_pos(:)),1);...
                vec;zeros(snum.idx(6,2)-length(vec),1)];
            %}
            % belief dynamics (note: this one could be written using
            % A_eq,b_eq, but it's error-prone and time-consuming. 
            % So I use this form.
            switch plan_mode
                case 'lin'
                    constrLinEq = @(s) this.getLinEqConstr_kf(s,snum,tar);
                case 'nl'
                    constrLinEq = @(s) this.getLinEqConstr(s,snum,tar);
            end
            
            %%% linear inequality constraints
            % bounds on states and input
            % this.w_lb <= u(1,:) <= this.w_ub;
            % this.a_lb <= u(2,:) <= this.a_ub;
            % this.v_lb <= z(4,:) <= this.v_ub;
            % [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];
            %{
            mat = diag([ones(snum.idx(2,2),1);zeros(snum.idx(6,2)-snum.idx(2,2),1)]); 
            idx = 3:4:(3+4*N); % since no constraint on z(3,:), set element to be 0
            for ii = idx
                mat(ii,ii) = 0;
            end
            A_ineq = [mat;-mat];
            vec1 = [repmat([fld.fld_cor(2);fld.fld_cor(4);0;this.v_ub],N+1,1);...
                repmat([this.w_ub;this.a_ub],N,1)];
            vec2 = -[repmat([fld.fld_cor(1);fld.fld_cor(3);0;this.v_lb],N+1,1);...
                repmat([this.w_lb;this.a_lb],N,1)];
            b_ineq = [vec1;vec2];
            %}
            
            constrLinIneq = @(s) this.getLinIneqConstr(s,snum,fld);
            
            penalty_coeff = cfg.initial_penalty_coeff; % Coefficient of l1 penalties 
            trust_box_size = cfg.initial_trust_box_size; % The trust region will be a box around the current iterate x.            
            hinge = @(x) sum(max(x,0));
            abssum = @(x) sum(abs(x));
            
            gam_iter = 1;
            %% loop 1: change alpha in \gamma modeling
            while(1)        
                fprintf('  [gamma loop] Robot.m, line %d.\n',MFileLineNr())
                fprintf('  gamma loop #%d\n',gam_iter)
                
                switch plan_mode
                    case 'lin'
                        objBel = @(s) this.getBelConstr_kf(s,snum,alp1,alp2,alp3); % belief update
                    case 'nl'
                        objBel = @(s) this.getBelConstr(s,snum,Kref,alp1,alp2,alp3); % belief update
                end
                
                %%% objective
                obj = @(s) this.getObj(s,snum,alp1,alp2,alp3);
                
                %% loop 2: penalty iteration
%                 ctol = 1e-2;
%                 mu = 0.1;
                pen_iter = 0;
                while(1)
                    pen_iter = pen_iter+1;
                    fprintf('    [Penalty loop] Robot.m, line %d\n', MFileLineNr())                                        
                    fprintf('    pentaly loop #%d\n',pen_iter);                    
                    
                    objNLineq = @(s) 0; % nonlinear inequality constraint here
                    objNLeq = @(s) [objKin(s);objBel(s)]; % nonlinear equality constraint here
%                     objNLeq = @(s) 0;
                    % trust region for approximating gamma. Trust region is the
                    % stepsize of change, i.e. abs(z-zref)
%                     tr = zeros(6,N); % bound of stepsize
%                     % z
%                     tr(1:2,:) = 3*ones(2,N);
%                     tr(3,:) = 3*ones(1,N);%pi/5*ones(1,N);
%                     tr(4,:) = 3*ones(1,N);
%                     % u
%                     tr(5:6,:) = 3*ones(2,N);
                    % x
                    
                    %% loop 3: trust region SQP
%                     thr_l = 1/4;
%                     thr_h = 1/2; %3/4;
%                     xtol = 0.1;
%                     infea_flag = false; % flagging whether cvx is infeasible
                    A_ineq = [];
                    b_ineq = [];
                    A_eq = [];
                    b_eq = [];
                    
                    % make initial solution feasible
                    [s,success] = this.find_closest_feasible_point(s,constrLinIneq,constrLinEq);
                    if (~success)
                        return;
                    end
                    
                    [s, trust_box_size, success, merit, model_merit, new_merit] = this.minimize_merit_function(s, fQ, q, obj, A_ineq, b_ineq, A_eq, b_eq, constrLinIneq, constrLinEq, objNLineq, objNLeq, hinge, abssum, cfg, penalty_coeff, trust_box_size, snum, fld);
                    % loop 3 ends
                    %{
                    while (1)
                        % robot state and control
                        cvx_begin sdp
                        variables z(4,N+1) u(2,N) x(2*this.gmm_num,N+1)
                        variable P(2,2,this.gmm_num,N+1) semidefinite % symmetric
                        % debug purpose
                        variable x_pred(2*this.gmm_num,N)
                        variable P_pred(2,2,this.gmm_num,N) semidefinite% symmetric
                        
                        % auxiliary variable
                        variable t(this.gmm_num*this.gmm_num,N+1)
                        variables slk_kin(4,N) slk_P(2,2,this.gmm_num,N) %slk_P_pred(2,2,this.gmm_num,N) 
                        expression t_unscaled(this.gmm_num*this.gmm_num,N+1)
                        expression obj
                        expression delta_x(2*this.gmm_num,N)
                        expression delta_P(2,2,this.gmm_num,N+1)
                        
                        % obj
                        delta_x = x-xref;
                        delta_P = P-Pref;
                        [grad,hess] = this.numerical_grad_hess(xref,Pref,full_hessian);
                        
                        obj = this.cmpObj(xref,Pref)+grad*[delta_x(:);delta_P(:)];%+...
%                             [delta_x(:);delta_P(:)]'*hess*[delta_x(:);delta_P(:)]/2;
                        obj = obj + mu*(sum(abs(slk_kin(:)))+...
                            sum(abs(slk_P(:)))); % add slack varaible %sum(abs(slk_P_pred(:)))+
                        
                        minimize(obj)
                        
                        % constraints
                        % initial value
                        z(:,1) == this.state;
                        x(:,1) == this.est_pos(:);
                        for jj = 1:this.gmm_num
                            triu(P(:,:,jj,1)) == triu(this.P{jj});
                        end
                        
                        % constraints on the go
                        for ii = 1:N
                            % linearize using previous result
                            % robot state
                            z(:,ii+1) == z(:,ii)+...
                                [z(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(z(3,ii)-zref(3,ii));
                                z(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(z(3,ii)-zref(3,ii));
                                u(:,ii)]*dt+slk_kin(:,ii);
                            [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];                           
                            
                            % target prediction
                            for jj = 1:this.gmm_num
                                %%%%% note: this part may need change later. In
                                %%%%% fact, linearziation should be wrt
                                %%%%% reference values.
                                A = del_f(x(2*jj-1:2*jj,ii));
                                if isempty (zref)
                                    C = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                                else
                                    C = del_h(xref(2*jj-1:2*jj,ii+1),zref(1:2,ii+1));
                                end
                                
                                % forward prediction
                                % mean
                                x_pred(2*jj-1:2*jj,ii) == f(x(2*jj-1:2*jj,ii));
                                % covariance
                                triu(P_pred(:,:,jj,ii)) == triu(A*P(:,:,jj,ii)*A'+Q);%+triu(slk_P_pred(:,:,jj,ii));
                                
                                % mean
                                %%%%% note: for now, I assume the mean is not
                                %%%%% affected by measurement in planning
                                x(2*jj-1:2*jj,ii+1) == x_pred(2*jj-1:2*jj,ii);
                                
                                % covariance
                                theta_bar = zeros(this.gmm_num,N+1);
                                for ll = 1:this.gmm_num
                                    tmp_vec = xref(2*ll-1:2*ll,:)-zref(1:2,:);
                                    theta_bar(ll,:) = atan2(tmp_vec(1,:),tmp_vec(2,:));
                                end
                                T = Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C;
                                expression tmp(this.gmm_num,2,2)
                                for ll = 1:this.gmm_num
                                    %%% note: gamma depends on ll, C depends
                                    %%% only on jj
                                    tmp(ll,1,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                                    tmp(ll,1,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(1,2,jj,ii),P_pred(2,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,2,jj,ii),P_pred_ref(2,2,jj,ii),alp1,alp2,alp3);
                                    %                             tmp(ll,2,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                    %                                 P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),theta_bar(ll,ii+1),...
                                    %                                 T(2,1),T(2,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii));
                                    tmp(ll,2,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(2,2,jj,ii),P_pred(1,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(2,2),T(2,1),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(2,2,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                                end
                                
                                triu(P(:,:,jj,ii+1)) == triu(squeeze(sum(tmp,1)))+triu(slk_P(2,2,jj,ii));
                            end
                        end
                        
                        this.w_lb <= u(1,:) <= this.w_ub;
                        this.a_lb <= u(2,:) <= this.a_ub;
                        this.v_lb <= z(4,:) <= this.v_ub;
                        
                         % trust region constraints
                         [-tr(1:4,ii) <= z(:,ii+1)-zref(:,ii+1) <= tr(1:4,ii)];
                         [-tr(5:6,ii) <= u(:,ii)-uref(:,ii) <= tr(5:6,ii)];
                        cvx_end
                        
                        % if infeasbile, use solution from last iteration                        
                        if strcmp(cvx_status,'Infeasible')
                            cprintf('Red',sprintf('Robot.m, line %d. Infeasibility issue\n',MFileLineNr()))
%                             error('inf')
                            infea_flag = true;
                            break
                        end
                        
                        %% determine whether to change the trust region
                        %
                        % compare the change in the objective function.
                        % true objective values
                        
                        % if P is nearly singluar, make it psd
                        for ii = 2:N+1
                            tmp = 0;
                            for ll = 1:this.gmm_num
                                P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/ 2;
                                mineigv = min(eig(P(:,:,ll,ii)));
                                if mineigv <= 0
                                    P(:,:,ll,ii) = P(:,:,ll,ii) + (abs(mineigv)+0.1)*eye(2);
                                end
                            end
                        end
                                                
                        act_prev = this.cmpMerit(zref,uref,zref,xref,Pref,P_pred_ref,Kref,mu);%this.cmpObj(xref,Pref); %this.cmpObj(xref,Pref);
                        act_cur = this.cmpMerit(z,u,zref,x,P,P_pred,Kref,mu);%this.cmpObj(x,P)
                        pred_prev = this.cmpMerit(zref,uref,zref,xref,Pref,P_pred_ref,Kref,mu); %this.cmpObj(xref,Pref)
                        pred_cur = cvx_optval;
                        
                        sprintf('Robot.m, line %d',MFileLineNr())
                        display('sqp objective value:')
                        display(cvx_optval)
                        display('actual objective value:')
                        display(act_cur)
                        % }
                        impv_ratio = (act_prev-act_cur)/(pred_prev-pred_cur);
                        display('improvement')
                        display(impv_ratio)
                        
                        % compare approximation with actual value
                        t = -5:0.1:5;
                        actval = zeros(length(t),1);
                        apprxval = zeros(length(t),1);
                        for pp = 1:length(t)
                            tt = t(pp);
                            ztmp = zref+tt;
                            utmp = uref+tt;
                            xtmp = xref+tt;
                            Ptmp = Pref+tt;
                            P_pred_tmp = P_pred_ref+tt;
                            
                            actval(pp) = this.cmpMerit(ztmp,utmp,zref,xtmp,Ptmp,P_pred_tmp,Kref,mu);
                            apprxval(pp) = this.cmpObj(xtmp,Ptmp)+grad*[xtmp(:)-xref(:);Ptmp(:)-Pref(:)];%+...
%                                 [xtmp(:)-xref(:);Ptmp(:)-Pref(:)]'*hess*[xtmp(:)-xref(:);Ptmp(:)-Pref(:)]/2;
                            apprxval(pp) = apprxval(pp) + mu*(sum(abs(slk_kin(:)))+...
                                sum(abs(slk_P(:))));
                        end
                        figure
                        hold on
                        plot(t,actval,'b')
                        plot(t,apprxval,'r')
                        
                        if pred_prev-pred_cur < -1e-5
                            display('approximate merit function got worse')
                            return
                        elseif impv_ratio < thr_h
                            % in this case, <Num Opt> book discusses two different
                            % cases, in one case, the new value is used. In
                            % another case, the variable values use the last
                            % iteration's. Here I just treat the varaible
                            % values unchanged
                            tr = tr*tr_dec;
                            cprintf('Green',sprintf('Robot.m, line %d. trust region shrinked\n',MFileLineNr()))
                            display(tr)
                            if all(abs(tr)<xtol)
                                break
                            end
                        elseif impv_ratio > thr_h
                            % determine if the optimal solution reaches the
                            % trust region boundary
                            tmp_z_dif = z(:,2:end)-zref(:,2:end);
                            tmp_u_dif = u-uref;
                            if max(abs([tmp_z_dif(:);tmp_u_dif(:)]-tr(:))) < 1e-4
                                % boundary reached
                                tr = tr*tr_inc;
                                cprintf('Blue',sprintf('Robot.m, line %d. trust region enlarged\n',MFileLineNr()))
                            else
                                cprintf('Cyan',sprintf('Robot.m, line %d. trust region unchanged\n',MFileLineNr()))
                            end
                            % variable values use the newly computed ones
                            zref = z;
                            uref = u;
                            xref = x;
                            Pref = P;
                            P_pred_ref = P_pred;
                            
                            %%% note sure if this is necessary, but I encounter cases
                            %%% where P_pred contains singular or non-psd terms. Can
                            %%% add the min eigvalue to make all terms psd. Haven't written code yet
                            
                            % compute Kref using Ricatti equation
                            for ii = 1:N
                                for jj = 1:this.gmm_num
                                    C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
                                    Kref(2*jj-1:2*jj,2*ii-1:2*ii) = P_pred(:,:,jj,ii)*C'/(C*P_pred(:,:,jj,ii)*C'+R);
                                end
                            end
                            break % break only when trust region is expanded
                        end
                    end % loop 3 ends
                    %}
                    if(hinge(objNLineq(s)) + abssum(objNLeq(s)) < cfg.cnt_tolerance || pen_iter >= cfg.max_penalty_iter ) %cfg.max_iter
                        
                        break;
                    end
                    trust_box_size = cfg.initial_trust_box_size;
                    penalty_coeff = cfg.merit_coeff_increase_ratio*penalty_coeff;
%                     if (any(abs(slk_kin(:))>ctol)... %||any(abs(slk_P_pred(:))>ctol)
%                             ||any(abs(slk_P(:))>ctol)) && ~infea_flag
%                         mu = mu*mu_inc;
%                     else
%                         break
%                     end
                end % loop 2 ends
                
                if ~success %infea_flag
                    % if CVS infeasible, directly reuse solution from
                    % previous step
                    break
                end
                
                x = this.convState(s,snum,'x');
%                 uref = this.convState(s,snum,'u');
                z = this.convState(s,snum,'z');
                
                % here we use the difference between the actual gamma and
                % gamma_aprx to decide the region.
                is_in_fov = zeros(this.gmm_num,N);
                gamma_exact = zeros(this.gmm_num,N);
%                 gamma_aprx = zeros(this.gmm_num,N);
                tmp_rbt = this;
                
                for ii = 1:N
                    for jj = 1:this.gmm_num
                        %%% this part can be revised when using probability
                        %%% of inFOV later
                        
                        tar_pos = x(2*jj-1:2*jj,ii+1); % use each gmm component mean as a possible target position
                        % actual inFOV
                        tmp_rbt.state = z(:,ii+1);
                        is_in_fov(jj,ii) = tmp_rbt.inFOV(tar_pos);
                        
                        % exact gamma
                        gamma_exact(jj,ii) = this.gam(z(1:2,ii+1),z(3,ii+1),...
                            tar_pos,alp1,alp2,alp3);
                        
                        % approximated inFOV
%                         gamma_aprx(jj,ii) = gam_aprx(z(1:2,ii+1),z(3,ii+1),...
%                             tar_pos,zref(1:2,ii+1),zref(3,ii+1),alp1,alp2,alp3);
                    end
                end
                
                %                     tmp_ratio = zeros(this.gmm_num,N);
                tmp_dif = zeros(this.gmm_num,N);
                for ii = 1:N
                    for jj = 1:this.gmm_num
%                         tmp_dif(jj,ii) = abs(gamma_aprx(jj,ii)-gamma_exact(jj,ii));
                        tmp_dif(jj,ii) = abs(is_in_fov(jj,ii)-gamma_exact(jj,ii));
                        %                             tmp_ratio(jj,ii) = abs((gamma_aprx(jj,ii)-gamma_exact(jj,ii))/max(gamma_exact(jj,ii),0.001));
                    end
                end
                
%                 tmp_z_dif = zeros(3,N+1);
%                 inc_flag = false(3,1); % if true, increase the trust region
                
                % update initial solution for next iteration
                s = this.updOptVar(s,snum,alp1,alp2,alp3);
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                if max(tmp_dif) <= cfg.gamma_tol %0.05
                    fprintf('  [gamma loop] Robot.m, line %d, robot state:\n', MFileLineNr())
                    fprintf('  Gamma error is within tolerance, break out of gamma loop\n')
                    break
                elseif gam_iter >= cfg.max_gam_iter 
                    fprintf('  [gamma loop] Robot.m, line %d, robot state:\n', MFileLineNr())
                    fprintf('  max gamma iteration reached. break out of gamma loop\n')
                    break
                else
                    cprintf('Magenta',sprintf('  [gamma loop] Robot.m, line %d.  gamma_exact is not close, go to another iteration.\n',MFileLineNr()))
%                     display(is_in_fov)
%                     display(gamma_exact)
                    alp1 = alp1*alp_inc;
                    alp2 = alp2*alp_inc;
                    alp3 = alp3*alp_inc;                    
                end            
                gam_iter = gam_iter+1;
            end % loop 1 ends
            
            uref = this.convState(s,snum,'u');
            
            if ~success %infea_flag
                uref = [uref(:,2:end),uref(:,end)];
            end

            % use actual dynamics to simulate
            zref = this.simState(uref,zref(:,1));
            optz = zref;
            optu = uref;
%             optz = zref;
%             optu = uref;
            
            % visualize the planned path
            %%% xref in this part needs change when infeasibility happens
            this.plotPlannedTraj(optz,xref,fld)
%             }
         end
        
        
        function [optz,optu] = cvxPlanner_sqp(this,fld,optz,optu) % cvxPlanner(this,fld,init_sol)
            % use the multi-layer approach similar to Sachin's work. Fix
            % the parameter for the sensor, solve path planning. Then
            % refine the parameter until close to reality. In each
            % iteration, a QP program is solved. The initial solution
            % comes from ngPlanner
            
            % planing in non-Gaussian (GMM) belief space
            N = this.mpc_hor;
            dt = this.dt;
            
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            alp_inc = 2; % increament paramter for alphatr_inc
            gam = @this.gam;
            gam_aprx = @this.gam_aprx;
           	p_aprx = @this.p_aprx;
            
            % trust region, penalty factor
            tr_inc = this.tr_inc;
            tr_dec = this.tr_dec;
            mu_inc = this.mu_inc;
            
            % set up simulation            
            if isempty(optz)
                prev_state = [];
            else
                prev_state = struct('optz',optz,'optu',optu);
            end
            
            init_sol = genInitState(this,fld,prev_state);
            
            zref = init_sol.z;
            uref = init_sol.u;
            xref = init_sol.x;
            Kref = init_sol.K;
            Pref = init_sol.P;
            P_pred_ref = init_sol.P_pred;
            full_hessian = false; 
            
            %% loop 1: change alpha in \gamma modeling
            while(1)
                %% loop 2: penalty iteration
                ctol = 1e-2;
                mu = 0.1;
                while(1)
                    % trust region for approximating gamma. Trust region is the
                    % stepsize of change, i.e. abs(z-zref)
                    tr = zeros(6,N); % bound of stepsize
                    % z
                    tr(1:2,:) = 3*ones(2,N);
                    tr(3,:) = 3*ones(1,N);%pi/5*ones(1,N);
                    tr(4,:) = 3*ones(1,N);
                    % u
                    tr(5:6,:) = 3*ones(2,N);
                    % x
                    
                    
                    
                    %                 for ii = 1:N
                    %                     for jj = 1:3
                    %                         tmp = 0.1*zref(jj,ii+1);
                    %                         if tmp == 0
                    %                             % the trust region of heading change is smaller than
                    %                             % the trust region for position
                    %                             if jj == 3
                    %                                 tmp = 0.1;
                    %                             else
                    %                                 tmp = 1;
                    %                             end
                    %                         end
                    %                         tr(jj,ii) = abs(tmp);
                    %                     end
                    %                 end
                    
                    %% loop 3: trust region SQP
                    thr_l = 1/4;
                    thr_h = 1/2; %3/4;
                    xtol = 0.1;
                    infea_flag = false; % flagging whether cvx is infeasible
                    while (1)
                        % robot state and control
                        cvx_begin sdp
                        variables z(4,N+1) u(2,N) x(2*this.gmm_num,N+1)
                        variable P(2,2,this.gmm_num,N+1) semidefinite % symmetric
                        % debug purpose
                        variable x_pred(2*this.gmm_num,N)
                        variable P_pred(2,2,this.gmm_num,N) semidefinite% symmetric
                        
                        % auxiliary variable
                        variable t(this.gmm_num*this.gmm_num,N+1)
                        variables slk_kin(4,N) slk_P(2,2,this.gmm_num,N) %slk_P_pred(2,2,this.gmm_num,N) 
                        expression t_unscaled(this.gmm_num*this.gmm_num,N+1)
                        expression obj
                        expression delta_x(2*this.gmm_num,N)
                        expression delta_P(2,2,this.gmm_num,N+1)
                        
                        % obj
                        delta_x = x-xref;
                        delta_P = P-Pref;
                        [grad,hess] = this.numerical_grad_hess(xref,Pref,full_hessian);
                        
                        obj = this.cmpObj(xref,Pref)+grad*[delta_x(:);delta_P(:)];%+...
%                             [delta_x(:);delta_P(:)]'*hess*[delta_x(:);delta_P(:)]/2;
                        obj = obj + mu*(sum(abs(slk_kin(:)))+...
                            sum(abs(slk_P(:)))); % add slack varaible %sum(abs(slk_P_pred(:)))+
                        
                        minimize(obj)
                        
                        % constraints
                        % initial value
                        z(:,1) == this.state;
                        x(:,1) == this.est_pos(:);
                        for jj = 1:this.gmm_num
                            triu(P(:,:,jj,1)) == triu(this.P{jj});
                        end
                        
                        % constraints on the go
                        for ii = 1:N
                            % linearize using previous result
                            % robot state
                            z(:,ii+1) == z(:,ii)+...
                                [z(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(z(3,ii)-zref(3,ii));
                                z(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(z(3,ii)-zref(3,ii));
                                u(:,ii)]*dt+slk_kin(:,ii);
                            [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];                           
                            
                            % target prediction
                            for jj = 1:this.gmm_num
                                %%%%% note: this part may need change later. In
                                %%%%% fact, linearziation should be wrt
                                %%%%% reference values.
                                A = del_f(x(2*jj-1:2*jj,ii));
                                if isempty (zref)
                                    C = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                                else
                                    C = del_h(xref(2*jj-1:2*jj,ii+1),zref(1:2,ii+1));
                                end
                                
                                % forward prediction
                                % mean
                                x_pred(2*jj-1:2*jj,ii) == f(x(2*jj-1:2*jj,ii));
                                % covariance
                                triu(P_pred(:,:,jj,ii)) == triu(A*P(:,:,jj,ii)*A'+Q);%+triu(slk_P_pred(:,:,jj,ii));
                                
                                % mean
                                %%%%% note: for now, I assume the mean is not
                                %%%%% affected by measurement in planning
                                x(2*jj-1:2*jj,ii+1) == x_pred(2*jj-1:2*jj,ii);
                                
                                % covariance
                                theta_bar = zeros(this.gmm_num,N+1);
                                for ll = 1:this.gmm_num
                                    tmp_vec = xref(2*ll-1:2*ll,:)-zref(1:2,:);
                                    theta_bar(ll,:) = atan2(tmp_vec(1,:),tmp_vec(2,:));
                                end
                                T = Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C;
                                expression tmp(this.gmm_num,2,2)
                                for ll = 1:this.gmm_num
                                    %%% note: gamma depends on ll, C depends
                                    %%% only on jj
                                    tmp(ll,1,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                                    tmp(ll,1,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(1,2,jj,ii),P_pred(2,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,2,jj,ii),P_pred_ref(2,2,jj,ii),alp1,alp2,alp3);
                                    %                             tmp(ll,2,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                    %                                 P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),theta_bar(ll,ii+1),...
                                    %                                 T(2,1),T(2,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii));
                                    tmp(ll,2,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                        P_pred(2,2,jj,ii),P_pred(1,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                        T(2,2),T(2,1),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(2,2,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                                end
                                
                                triu(P(:,:,jj,ii+1)) == triu(squeeze(sum(tmp,1)))+triu(slk_P(2,2,jj,ii));
                            end
                        end
                        
                        this.w_lb <= u(1,:) <= this.w_ub;
                        this.a_lb <= u(2,:) <= this.a_ub;
                        this.v_lb <= z(4,:) <= this.v_ub;
                        
                         % trust region constraints
                         [-tr(1:4,ii) <= z(:,ii+1)-zref(:,ii+1) <= tr(1:4,ii)];
                         [-tr(5:6,ii) <= u(:,ii)-uref(:,ii) <= tr(5:6,ii)];
                        cvx_end
                        
                        % if infeasbile, use solution from last iteration                        
                        if strcmp(cvx_status,'Infeasible')
                            cprintf('Red',sprintf('Robot.m, line %d. Infeasibility issue\n',MFileLineNr()))
%                             error('inf')
                            infea_flag = true;
                            break
                        end
                        
                        %% determine whether to change the trust region
                        %
                        % compare the change in the objective function.
                        % true objective values
                        
                        % if P is nearly singluar, make it psd
                        for ii = 2:N+1
                            tmp = 0;
                            for ll = 1:this.gmm_num
                                P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/ 2;
                                mineigv = min(eig(P(:,:,ll,ii)));
                                if mineigv <= 0
                                    P(:,:,ll,ii) = P(:,:,ll,ii) + (abs(mineigv)+0.1)*eye(2);
                                end
                            end
                        end
                                                
                        act_prev = this.cmpMerit(zref,uref,zref,xref,Pref,P_pred_ref,Kref,mu);%this.cmpObj(xref,Pref); %this.cmpObj(xref,Pref);
                        act_cur = this.cmpMerit(z,u,zref,x,P,P_pred,Kref,mu);%this.cmpObj(x,P)
                        pred_prev = this.cmpMerit(zref,uref,zref,xref,Pref,P_pred_ref,Kref,mu); %this.cmpObj(xref,Pref)
                        pred_cur = cvx_optval;
                        
                        sprintf('Robot.m, line %d',MFileLineNr())
                        display('sqp objective value:')
                        display(cvx_optval)
                        display('actual objective value:')
                        display(act_cur)
                        %}
                        impv_ratio = (act_prev-act_cur)/(pred_prev-pred_cur);
                        display('improvement')
                        display(impv_ratio)
                        
                        % compare approximation with actual value
                        t = -5:0.1:5;
                        actval = zeros(length(t),1);
                        apprxval = zeros(length(t),1);
                        for pp = 1:length(t)
                            tt = t(pp);
                            ztmp = zref+tt;
                            utmp = uref+tt;
                            xtmp = xref+tt;
                            Ptmp = Pref+tt;
                            P_pred_tmp = P_pred_ref+tt;
                            
                            actval(pp) = this.cmpMerit(ztmp,utmp,zref,xtmp,Ptmp,P_pred_tmp,Kref,mu);
                            apprxval(pp) = this.cmpObj(xtmp,Ptmp)+grad*[xtmp(:)-xref(:);Ptmp(:)-Pref(:)];%+...
%                                 [xtmp(:)-xref(:);Ptmp(:)-Pref(:)]'*hess*[xtmp(:)-xref(:);Ptmp(:)-Pref(:)]/2;
                            apprxval(pp) = apprxval(pp) + mu*(sum(abs(slk_kin(:)))+...
                                sum(abs(slk_P(:))));
                        end
                        figure
                        hold on
                        plot(t,actval,'b')
                        plot(t,apprxval,'r')
                        
                        if pred_prev-pred_cur < -1e-5
                            display('approximate merit function got worse')
                            return
                        elseif impv_ratio < thr_h
                            % in this case, <Num Opt> book discusses two different
                            % cases, in one case, the new value is used. In
                            % another case, the variable values use the last
                            % iteration's. Here I just treat the varaible
                            % values unchanged
                            tr = tr*tr_dec;
                            cprintf('Green',sprintf('Robot.m, line %d. trust region shrinked\n',MFileLineNr()))
                            display(tr)
                            if all(abs(tr)<xtol)
                                break
                            end
                        elseif impv_ratio > thr_h
                            % determine if the optimal solution reaches the
                            % trust region boundary
                            tmp_z_dif = z(:,2:end)-zref(:,2:end);
                            tmp_u_dif = u-uref;
                            if max(abs([tmp_z_dif(:);tmp_u_dif(:)]-tr(:))) < 1e-4
                                % boundary reached
                                tr = tr*tr_inc;
                                cprintf('Blue',sprintf('Robot.m, line %d. trust region enlarged\n',MFileLineNr()))
                            else
                                cprintf('Cyan',sprintf('Robot.m, line %d. trust region unchanged\n',MFileLineNr()))
                            end
                            % variable values use the newly computed ones
                            zref = z;
                            uref = u;
                            xref = x;
                            Pref = P;
                            P_pred_ref = P_pred;
                            
                            %%% note sure if this is necessary, but I encounter cases
                            %%% where P_pred contains singular or non-psd terms. Can
                            %%% add the min eigvalue to make all terms psd. Haven't written code yet
                            
                            % compute Kref using Ricatti equation
                            for ii = 1:N
                                for jj = 1:this.gmm_num
                                    C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
                                    Kref(2*jj-1:2*jj,2*ii-1:2*ii) = P_pred(:,:,jj,ii)*C'/(C*P_pred(:,:,jj,ii)*C'+R);
                                end
                            end
                            break % break only when trust region is expanded
                        end
                    end % loop 3 ends
                    
                    if (any(abs(slk_kin(:))>ctol)... %||any(abs(slk_P_pred(:))>ctol)
                            ||any(abs(slk_P(:))>ctol)) && ~infea_flag
                        mu = mu*mu_inc;
                    else
                        break
                    end
                end % loop 2 ends
                
                if infea_flag
                    % if CVS infeasible, directly reuse solution from
                    % previous step
                    break
                end
                
                % here we use the difference between the actual gamma and
                % gamma_aprx to decide the region.
                is_in_fov = zeros(this.gmm_num,N);
                gamma_exact = zeros(this.gmm_num,N);
                gamma_aprx = zeros(this.gmm_num,N);
                tmp_rbt = this;
                for ii = 1:N
                    for jj = 1:this.gmm_num
                        %%% this part can be revised when using probability
                        %%% of inFOV later
                        
                        tar_pos = x(2*jj-1:2*jj,ii+1); % use each gmm component mean as a possible target position
                        % actual inFOV
                        tmp_rbt.state = z(:,ii+1);
                        is_in_fov(jj,ii) = tmp_rbt.inFOV(tar_pos);
                        
                        % exact gamma
                        gamma_exact(jj,ii) = gam(z(1:2,ii+1),z(3,ii+1),...
                            tar_pos,alp1,alp2,alp3);
                        
                        % approximated inFOV
                        gamma_aprx(jj,ii) = gam_aprx(z(1:2,ii+1),z(3,ii+1),...
                            tar_pos,zref(1:2,ii+1),zref(3,ii+1),alp1,alp2,alp3);
                    end
                end
                
                %                     tmp_ratio = zeros(this.gmm_num,N);
                tmp_dif = zeros(this.gmm_num,N);
                for ii = 1:N
                    for jj = 1:this.gmm_num
%                         tmp_dif(jj,ii) = abs(gamma_aprx(jj,ii)-gamma_exact(jj,ii));
                        tmp_dif(jj,ii) = abs(is_in_fov(jj,ii)-gamma_exact(jj,ii));
                        %                             tmp_ratio(jj,ii) = abs((gamma_aprx(jj,ii)-gamma_exact(jj,ii))/max(gamma_exact(jj,ii),0.001));
                    end
                end
                
%                 tmp_z_dif = zeros(3,N+1);
%                 inc_flag = false(3,1); % if true, increase the trust region
                
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                if max(tmp_dif) <= 0.05
                    break
                else
                    cprintf('Magenta',sprintf('Robot.m, line %d.  gamma_exact is not close',MFileLineNr()))
                    display(is_in_fov)
                    display(gamma_exact)
                    alp1 = alp1*alp_inc;
                    alp2 = alp2*alp_inc;
                    alp3 = alp3*alp_inc;
                end                
            end % loop 1 ends
            
            if infea_flag
                uref = [uref(:,2:end),uref(:,end)];
            end

            % use actual dynamics to simulate
            zref = this.simState(uref);
            optz = zref;
            optu = uref;
%             optz = zref;
%             optu = uref;
            
            % visualize the planned path
            %%% xref in this part needs change when infeasibility happens
            this.plotPlannedTraj(optz,xref,fld)
%             }
        end
        
        function [optz,optu] = cvxPlanner(this,fld,optz,optu) % cvxPlanner(this,fld,init_sol)
            % use the multi-layer approach similar to Sachin's work. Fix
            % the parameter for the sensor, solve path planning. Then
            % refine the parameter until close to reality. In each
            % iteration, a convex program is solved. The initial solution
            % comes from ngPlanner
            % note this formulation is problematic since the
            % log-determinant term in obj makes the obj non-cvx. Convexity
            % holds iff we treat covariance as a constant, using P_ref.
            
            % planing in non-Gaussian (GMM) belief space
            N = this.mpc_hor;
            dt = this.dt;
            
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            alp_inc = 2; % increament paramter for alphatr_inc
            gam = @this.gam;
            gam_aprx = @this.gam_aprx;
           	p_aprx = @this.p_aprx;
            
            % trust region
            tr_inc = this.tr_inc;
            tr_dec = this.tr_dec;
            
            % set up simulation            
            if isempty(optz)
                prev_state = [];
            else
                prev_state = struct('optz',optz,'optu',optu);
            end
            
            init_sol = genInitState(this,fld,prev_state);
            
            zref = init_sol.z;
            uref = init_sol.u;
            xref = init_sol.x;
            Kref = init_sol.K;
            Pref = init_sol.P;
            P_pred_ref = init_sol.P_pred;

            % outer loop: change alpha in \gamma modeling
            while(1)
                % inner loop: change the trust region
                
                % trust region for approximating gamma. Trust region is the
                % stepsize of change, i.e. abs(z-zref)
                tr = zeros(3,N); % bound of stepsize
                tr(1:2,:) = 3*ones(2,N);
                tr(3,:) = pi/5*ones(1,N);
%                 for ii = 1:N
%                     for jj = 1:3
%                         tmp = 0.1*zref(jj,ii+1);
%                         if tmp == 0
%                             % the trust region of heading change is smaller than
%                             % the trust region for position
%                             if jj == 3
%                                 tmp = 0.1;
%                             else
%                                 tmp = 1;
%                             end
%                         end
%                         tr(jj,ii) = abs(tmp);
%                     end
%                 end
                
                while (1)
                    % robot state and control
                    cvx_begin sdp
                    variables z(4,N+1) u(2,N) x(2*this.gmm_num,N+1)
                    variable P(2,2,this.gmm_num,N+1) symmetric
                    % debug purpose
                    variable x_pred(2*this.gmm_num,N)
                    variable P_pred(2,2,this.gmm_num,N) symmetric
                    
                    % auxiliary variable
                    variable t(this.gmm_num*this.gmm_num,N+1)
                    expression t_unscaled(this.gmm_num*this.gmm_num,N+1)
                    expression obj
                                       
                    % epigraph for t
                    t(:) >= 0;
                    % LMI
                    for ii = 1:N
                        for jj = 1:this.gmm_num
                            for ll = 1:this.gmm_num
                                scale_factor = norm(xref(2*jj-1:2*jj,ii+1)-xref(2*ll-1:2*ll,ii+1));
                                if scale_factor == 0
                                    scale_factor = 1;
                                end
                                [P(:,:,ll,ii+1) (x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))/scale_factor;
                                    (x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))'/scale_factor t(this.gmm_num*(jj-1)+ll,ii+1)]>=0;
                                % because of the scale_factor, t should be
                                % scaled back to compute the obj. 
                                t_unscaled(this.gmm_num*(jj-1)+ll,ii+1) = t(this.gmm_num*(jj-1)+ll,ii+1)*(scale_factor)^2;
                            end
                        end
                    end
                    
                    % obj 
                    obj = 0;
                    for ii = 1:N
                        for jj = 1:this.gmm_num
                            % compute the -log(2pi*|P|)
                            tmp = 0;
%                             for ll = 1:this.gmm_num
%                                 tmp = tmp-log_det(P(:,:,ll,ii+1))-log(2*pi);
%                             end
                            % compute the whole objective
%                             obj = obj+this.wt(jj)*(sum(t_unscaled(this.gmm_num*(jj-1)+1:this.gmm_num*jj,ii+1))/2+tmp);                            
%                             obj = obj+this.wt(jj)*sum(t(this.gmm_num*(jj-1)+1:this.gmm_num*jj,ii+1));          
                            obj = obj+this.wt(jj)*sum(t_unscaled(this.gmm_num*(jj-1)+1:this.gmm_num*jj,ii+1))/2;     
                        end
                    end                    
                    
                    % obj
                    minimize(obj) %sum(t(:))
                    
                    % constraints
                    % initial value
                    z(:,1) == this.state;
                    x(:,1) == this.est_pos(:);
                    for jj = 1:this.gmm_num
                        triu(P(:,:,jj,1)) == triu(this.P{jj});
                    end
                    
                    % constraints on the go
                    for ii = 1:N
                        % linearize using previous result
                        % robot state
                        z(:,ii+1) == z(:,ii)+...
                            [z(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(z(3,ii)-zref(3,ii));
                            z(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(z(3,ii)-zref(3,ii));
                            u(:,ii)]*dt;
                        [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];
                        
                        % trust region constraints
                        [-tr(:,ii) <= z(1:3,ii+1)-zref(1:3,ii+1) <= tr(:,ii)];
                        
                        % target prediction
                        for jj = 1:this.gmm_num
                            A = del_f(x(2*jj-1:2*jj,ii));
                            if isempty (zref)
                                C = del_h(x(2*jj-1:2*jj,ii),z(1:2,ii));
                            else
                                C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
                            end
                            
                            % forward prediction
                            % mean
                            x_pred(2*jj-1:2*jj,ii) == f(x(2*jj-1:2*jj,ii));
                            % covariance
                            triu(P_pred(:,:,jj,ii)) == triu(A*P(:,:,jj,ii)*A'+Q);
                            
                            % mean
                            %%%%% note: for now, I assume the mean is not
                            %%%%% affected by measurement in planning
                            x(2*jj-1:2*jj,ii+1) == x_pred(2*jj-1:2*jj,ii);
                            
                            % covariance
                            theta_bar = zeros(this.gmm_num,N+1);
                            for ll = 1:this.gmm_num
                                tmp_vec = xref(2*ll-1:2*ll,:)-zref(1:2,:);
                                theta_bar(ll,:) = atan2(tmp_vec(1,:),tmp_vec(2,:));
                            end
                            T = Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C;
                            expression tmp(this.gmm_num,2,2)
                            for ll = 1:this.gmm_num
                                %%% note: gamma depends on ll, C depends
                                %%% only on jj
                                tmp(ll,1,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                    P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                    T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                                tmp(ll,1,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                    P_pred(1,2,jj,ii),P_pred(2,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                    T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,2,jj,ii),P_pred_ref(2,2,jj,ii),alp1,alp2,alp3);
                                %                             tmp(ll,2,1) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                %                                 P_pred(1,1,jj,ii),P_pred(2,1,jj,ii),xref(2*ll-1:2*ll,ii+1),theta_bar(ll,ii+1),...
                                %                                 T(2,1),T(2,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,jj,ii),P_pred_ref(2,1,jj,ii));
                                tmp(ll,2,2) = this.wt(ll)*p_aprx(z(1:2,ii+1),z(3,ii+1),...
                                    P_pred(2,2,jj,ii),P_pred(1,2,jj,ii),xref(2*ll-1:2*ll,ii+1),...
                                    T(2,2),T(2,1),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(2,2,jj,ii),P_pred_ref(2,1,jj,ii),alp1,alp2,alp3);
                            end
                            
                            triu(P(:,:,jj,ii+1)) == triu(squeeze(sum(tmp,1)));                            
                        end
                    end
                    this.w_lb <= u(1,:) <= this.w_ub;
                    this.a_lb <= u(2,:) <= this.a_ub;
                    this.v_lb <= z(4,:) <= this.v_ub;
                    
                    cvx_end
                    
                    % if infeasbile
                    if strcmp(cvx_status,'Infeasible')
                        sprintf('Robot.m, line %d',MFileLineNr())
                        error('inf')
                    end
                    
                    %% determine whether to change the trust region
                    %
                    
                    % compare the change in the objective function.
                    % true objective values
                    obj_prev = this.cmpObj(xref,Pref);
                    obj_cur = this.cmpObj(x,P);
                    % approximate objective values
                    obj_aprx_prev = this.cmpObjAprx(xref,Pref);
                    obj_aprx_cur = cvx_optval;
                    
                    % this snippet is used to test whether cvx obj is
                    % correctly formulated. cmpObjAprx gives the correct
                    % result. erros in cvx obj has been fixed by using
                    % cmpObjAprx. leave this code for record.
                    %{
                    optcvx = cvx_optval;
                    for ii = 1:N
                        for jj = 1:this.gmm_num
                            for ll = 1:this.gmm_num
                                optcvx = optcvx+this.wt(jj)*(log(det(P(:,:,ll,ii+1)))/2+log(2*pi));
                            end
                        end
                    end
                    %}
                    
                    sprintf('Robot.m, line %d',MFileLineNr())
%                     display(optcvx)
                    display('approximate objective value:/n')
                    display(this.cmpObjAprx(x,P))                                       
                    
                    % here we use the difference between the actual gamma and
                    % gamma_aprx to decide the region.
                    is_in_fov = zeros(this.gmm_num,N);
                    gamma_exact = zeros(this.gmm_num,N);
                    gamma_aprx = zeros(this.gmm_num,N);
                    tmp_rbt = this;
                    for ii = 1:N
                        for jj = 1:this.gmm_num
                            tar_pos = x(2*jj-1:2*jj,ii+1); % use each gmm component mean as a possible target position
                            % actual inFOV
                            tmp_rbt.state = z(:,ii+1);
                            is_in_fov(jj,ii) = tmp_rbt.inFOV(tar_pos);
                            
                            % exact gamma
                            gamma_exact(jj,ii) = gam(z(1:2,ii+1),z(3,ii+1),...
                                tar_pos,alp1,alp2,alp3);
                            
                            % approximated inFOV
                            gamma_aprx(jj,ii) = gam_aprx(z(1:2,ii+1),z(3,ii+1),...
                                tar_pos,zref(1:2,ii+1),zref(3,ii+1),alp1,alp2,alp3);
                        end
                    end
                    
%                     tmp_ratio = zeros(this.gmm_num,N);
                    tmp_dif = zeros(this.gmm_num,N);
                    for ii = 1:N
                        for jj = 1:this.gmm_num
                            tmp_dif(jj,ii) = abs(gamma_aprx(jj,ii)-gamma_exact(jj,ii));
%                             tmp_ratio(jj,ii) = abs((gamma_aprx(jj,ii)-gamma_exact(jj,ii))/max(gamma_exact(jj,ii),0.001));
                        end
                    end
                    
                    tmp_z_dif = zeros(3,N+1);
                    inc_flag = false(3,1); % if true, increase the trust region
                    
                    % if the difference between exact and estimated gamma
                    % is small, either increase the trust region or jump
                    % out the inner while loop
                    % if max(tmp_ratio) <= 0.1
                    if max(tmp_dif) <= 0.2
                        % determine if the optimal solution reaches the
                        % trust region boundary
                        for ii = 1:N+1
                            tmp_z_dif(:,ii) = z(1:3,ii)-zref(1:3,ii);
                        end
                        
                        for jj = 1:3
                            if max(abs(abs(tmp_z_dif(jj,:))-tr(jj))) <= 0.001
                                inc_flag(jj) = true;
                                break
                            end
                        end
                        
                        %                         for jj = 1:3
                        %                             if inc_flag(jj)
                        %                                 tr(jj,:) = tr(jj,:)*tr_inc;
                        %                             end
                        %                         end
                        
                        if sum(inc_flag) > 0
                            tr = tr*tr_inc;
                            cprintf('Red',sprintf('Robot.m, line %d. trust region enlarged\n',MFileLineNr()))
                        else
                            % if no increment is needed, break from the
                            % inner loop
                            cprintf('Cyan',sprintf('Robot.m, line %d. trust region unchanged\n',MFileLineNr()))
                            break
                        end
                        display(tr)
                        
                        % shrink trust region
                    else
                        tr = tr*tr_dec;
                        cprintf('Green',sprintf('Robot.m, line %d. trust region shrinked\n',MFileLineNr()))
                        display(tr)
                    end
                    
                end
                
                % assign value for next iteration
                zref = z;
                uref = u;
                xref = x;
                P_pred_ref = P_pred;
                
                %%% note sure if this is necessary, but I encounter cases
                %%% where P_pred contains singular or non-psd terms. Can
                %%% add the min eigvalue to make all terms psd. Haven't written code yet                
                
                
                % coompute Kref using Ricatti equation
                for ii = 1:N
                    for jj = 1:this.gmm_num
                        C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
                        Kref(2*jj-1:2*jj,2*ii-1:2*ii) = P_pred(:,:,jj,ii)*C'/(C*P_pred(:,:,jj,ii)*C'+R);
                    end
                end
                
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                dif = max(abs(is_in_fov-gamma_exact));
                if dif < 0.2
                    %                     display ('Robot.m line 829')
                    %                     display(P)
                    break
                else
                    cprintf('Blue',sprintf('Robot.m, line %d.  gamma_exact is not close',MFileLineNr()))
                    display(is_in_fov)
                    display(gamma_exact)
                end
                
                alp1 = alp1*alp_inc;
                alp2 = alp2*alp_inc;
                alp3 = alp3*alp_inc;
                
%                 display ('Robot.m line 522')
%                 display (alp1)
%                 display (alp2)
%                 display (alp3)
            end
            
            optz = zref;
            optu = uref;
            
            % visualize the planned path
            this.plotPlannedTraj(optz,xref,fld)
%             }
        end
        
        function [optz,optu] = cvxPlanner_kf(this,fld,optz,optu) % cvxPlanner(this,fld,init_sol)
            % this function is used with KF. The purpose is to test the
            % planning algorithm when using KF (a simpler case than gmm).
            
            % planing in non-Gaussian (GMM) belief space
            N = this.mpc_hor;
            dt = this.dt;
            
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            alp_inc = 2; % increament paramter for alpha
            gam = @this.gam;
            gam_aprx = @this.gam_aprx;
           	p_aprx = @this.p_aprx;
            
            % trust region
            tr_inc = this.tr_inc;
            tr_dec = this.tr_dec;
            
            % set up simulation             
            if isempty(optz)
                prev_state = [];
            else
                prev_state = struct('optz',optz,'optu',optu);
            end
            
            init_sol = genInitState_kf(this,fld,prev_state);
            
            zref = init_sol.z;
            uref = init_sol.u;
            xref = init_sol.x;
            Kref = init_sol.K;
            P_pred_ref = init_sol.P_pred;
            
            % outer loop: change alpha in \gamma modeling
            while (1)
                % inner loop: change the trust region
                
                % trust region for approximating gamma. Trust region is the
                % stepsize of change, i.e. abs(z-zref)
                tr = zeros(3,N); % bound of stepsize
                for ii = 1:N
                    for jj = 1:3
                        tmp = 0.1*zref(jj,ii+1);
                        if tmp == 0
                            % the trust region of heading change is smaller than
                            % the trust region for position
                            if jj == 3                                
                                tmp = 0.1;
                            else
                                tmp = 1;
                            end
                        end
                        tr(jj,ii) = abs(tmp);                      
                    end
                end
                
                while(1)
                    %% solving sequential cvx
                    % robot state and control
                    cvx_begin sdp
                    variables z(4,N+1) u(2,N) x(2,N+1)
                    variable P(2,2,N+1) semidefinite
                    % debug purpose
                    variable x_pred(2,N)
                    variable P_pred(2,2,N) semidefinite
                    
                    % auxiliary variable
                    expression t(N+1,1)
                    % obj
                    for ii = 1:N+1
                        t(ii) = trace(P(:,:,ii));
                    end
                    
                    minimize sum(t)
                    
                    % initial value
                    z(:,1) == this.state;
                    x(:,1) == this.est_pos(:);
                    triu(P(:,:,1)) == triu(this.P);
                    
                    % constraints on the go
                    for ii = 1:N
                        % linearize using previous result
                        % robot state
                        z(:,ii+1) == z(:,ii)+...
                            [z(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(z(3,ii)-zref(3,ii));
                            z(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(z(3,ii)-zref(3,ii));
                            u(:,ii)]*dt;
                        [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];
                        
                        % trust region constraints
                        [-tr(:,ii) <= z(1:3,ii+1)-zref(1:3,ii+1) <= tr(:,ii)];
                        
                        % target prediction
                        %                     for jj = 1:this.gmm_num
                        A = del_f(x(:,ii));
                        if isempty (zref)
                            C = del_h(x(:,ii),z(1:2,ii));
                        else
                            C = del_h(xref(:,ii),zref(1:2,ii));
                        end
                        
                        % forward prediction
                        % mean
                        x_pred(:,ii) == f(x(:,ii));
                        % covariance
                        triu(P_pred(:,:,ii)) == triu(A*P(:,:,ii)*A'+Q);
                        
                        % mean
                        %%%%% note: for now, I assume the mean is not
                        %%%%% affected by measurement in planning
                        x(:,ii+1) == x_pred(:,ii);
                        
                        % covariance
                        T = Kref(:,2*ii-1:2*ii)*C;
                        expression tmp(2,2)
                        tmp(1,1) = p_aprx(z(1:2,ii+1),z(3,ii+1),...
                            P_pred(1,1,ii),P_pred(2,1,ii),xref(:,ii+1),...
                            T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,1,ii),P_pred_ref(2,1,ii),alp1,alp2,alp3);
                        tmp(1,2) = p_aprx(z(1:2,ii+1),z(3,ii+1),...
                            P_pred(1,2,ii),P_pred(2,2,ii),xref(:,ii+1),...
                            T(1,1),T(1,2),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(1,2,ii),P_pred_ref(2,2,ii),alp1,alp2,alp3);
                        tmp(2,2) = p_aprx(z(1:2,ii+1),z(3,ii+1),...
                            P_pred(2,2,ii),P_pred(1,2,ii),xref(:,ii+1),...
                            T(2,2),T(2,1),zref(1:2,ii+1),zref(3,ii+1),P_pred_ref(2,2,ii),P_pred_ref(2,1,ii),alp1,alp2,alp3);
                        
                        triu(P(:,:,ii+1)) == triu(tmp);
                    end
                    this.w_lb <= u(1,:) <= this.w_ub;
                    this.a_lb <= u(2,:) <= this.a_ub;
                    this.v_lb <= z(4,:) <= this.v_ub;
                    
                    cvx_end
                   
                    % if infeasbile
                    if strcmp(cvx_status,'Infeasible')
                        sptrinf('Robot.m, line %d',MFileLineNr())
                        error('inf')
                    end
                    
                    %% determine whether to change the trust region 
                    % here we use the difference between the actual gamma and
                    % gamma_aprx to decide the region.
                    is_in_fov = zeros(N,1);
                    gamma_exact = zeros(N,1);
                    gamma_aprx = zeros(N,1);
                    tmp_rbt = this;
                    for ii = 1:N
                        tar_pos = x(:,ii+1); % use each gmm component mean as a possible target position
                        % actual inFOV
                        tmp_rbt.state = z(:,ii+1);
                        is_in_fov(ii) = tmp_rbt.inFOV(tar_pos);
                        
                        % exact gamma
                        gamma_exact(ii) = gam(z(1:2,ii+1),z(3,ii+1),...
                            tar_pos,alp1,alp2,alp3);
                        
                        % approximated inFOV
                        gamma_aprx(ii) = gam_aprx(z(1:2,ii+1),z(3,ii+1),...
                            tar_pos,zref(1:2,ii+1),zref(3,ii+1),alp1,alp2,alp3);
                        
                    end
                    
                    % ratio = (gamma_aprx-gamma_exact)/gamma_exact
                    tmp_ratio = zeros(N,1);
                    for ii = 1:N
                        tmp_ratio(ii) = abs((gamma_aprx(ii)-gamma_exact(ii))/max(gamma_exact(ii),0.001));
                    end
                    
                    tmp_z_dif = zeros(3,N+1);
                    inc_flag = false(3,1); % if true, increase the trust region                   
                    
                    % if the difference between exact and estimated gamma
                    % is small, either increase the trust region or jump
                    % out the inner while loop
                    if max(tmp_ratio) <= 0.1                        
                        % determine if the optimal solution reaches the 
                        % trust region boundary
                        for ii = 1:N+1
                            tmp_z_dif(:,ii) = z(1:3,ii)-zref(1:3,ii);                            
                        end
                        
                        for jj = 1:3
                           if max(abs(abs(tmp_z_dif(jj,:))-tr(jj))) <= 0.001
                                inc_flag(jj) = true;
                                break
                           end                                                       
                        end
                        
%                         for jj = 1:3
%                             if inc_flag(jj)
%                                 tr(jj,:) = tr(jj,:)*tr_inc;
%                             end
%                         end
                        if sum(inc_flag) > 0
                            tr = tr*tr_inc;
                            cprintf('Red',sprintf('Robot.m, line %d. trust region enlarged\n',MFileLineNr()))
                        else
                            % if no increment is needed, break from the
                            % inner loop
                            cprintf('Cyan',sprintf('Robot.m, line %d. trust region unchanged\n',MFileLineNr()))
                            break
                        end                        
                        display(tr)
                        
                    % shrink trust region
                    else
                        tr = tr*tr_dec;
                        cprintf('Green',sprintf('Robot.m, line %d. trust region shrinked\n',MFileLineNr()))
                        display(tr)
                    end
                    
                end
                
                % assign value for next iteration
                zref = z;
                uref = u;
                xref = x;                
                P_pred_ref = P_pred;
                Pref = P;
                
                % compute Kref using Ricatti equation 
                for ii = 1:N
                    C = del_h(xref(:,ii),zref(1:2,ii));
                    Kref(:,2*ii-1:2*ii) = P_pred(:,:,ii)*C'/(C*P_pred(:,:,ii)*C'+R);
                end
                
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                dif = max(abs(is_in_fov-gamma_exact));
                if dif < 0.2
%                     display ('Robot.m line 829')
%                     display(P)
                    break
                else
                    cprintf('Blue',sprintf('Robot.m, line %d.  gamma_exact is not close',MFileLineNr()))
                    display(is_in_fov)
                    display(gamma_exact)
                end                               
                
                alp1 = alp1*alp_inc;
                alp2 = alp2*alp_inc;
                alp3 = alp3*alp_inc;
                
%                 display ('Robot.m line 838')
%                 display (alp1)
%                 display (alp2)
%                 display (alp3)
            end
            
            optz = zref;
            optu = uref;
            
            % visualize the planned path
            this.plotPlannedTraj(optz,xref,fld)
%             }
        end
     
        function [optz,optu] = ngPlanner(this,fld,optz,optu)
            % use the multi-layer approach similar to Sachin's work. Fix
            % the parameter for the sensor, solve path planning. Then
            % refine the parameter until close to reality.
            
            % planing in non-Gaussian (GMM) belief space
            N = this.mpc_hor;
            dt = this.dt;
            
            % target
            tar = fld.target;
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            % sensor
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            alp1 = this.alp1;
            alp2 = this.alp2;
            gam = this.gam;
            gam_aprx = this.gam_aprx;
            gam_den = this.gam_den;
            
            % the parameter for the sensing boundary approximation
            alp = 1;
            alp_inc = 2; % increament paramter for alpha
            
            %%% extrapolate last-step solution as the seed solution for
            %%% current iteration. If this strategy does not work, then use
            %%% sequential cvx programming for all remaining planning work.
            %%%             
            % initial solution for approximation
            zref = [];
            x_ol_pred = zeros(2*this.gmm_num,N+1);
            x_ol_pred(:,1) = this.est_pos(:);
            for ii = 1:N
                x_ol_pred(:,ii+1) = f(x_ol_pred(:,ii));
            end
%             
%             theta_bar = zeros(this.gmm_num,1);
%             for jj = 1:this.gmm_num           
%                 tmp_vec = x_ol_pred(2*jj-1:2*jj,1)-this.state(1:2);
%                 theta_bar(jj) = atan2(tmp_vec(2),tmp_vec(1));
%             end
            if isempty(optz)
                prev_state = [];
            else
                prev_state = struct('optz',optz,'optu',optu);
            end
            
            init_sol = genInitState(this,fld,prev_state);
            % open-loop prediction of target
            x_olp = init_sol.x_olp;            
%             theta_bar = init_sol.theta_bar;
            
            % set up
            % robot state and control
            z = sdpvar(4,N+1,'full'); % robot state
            u = sdpvar(2,N,'full'); % robot control
            % estimation
            x = sdpvar(2*this.gmm_num,N+1,'full'); % target mean
            P = sdpvar(2,2,this.gmm_num,N+1,'full');
%             P = cell(this.gmm_num,N+1);
%             for ii = 1:N+1
%                 for jj = 1:this.gmm_num
%                     P{jj,ii} = sdpvar(2,2,'full'); % a set of 2-by-2 symmetric matrices
%                 end
%             end
            
            x_pred = sdpvar(2*this.gmm_num,N,'full'); % target mean
            P_pred = sdpvar(2,2,this.gmm_num,N,'full');
%             P_pred = cell(this.gmm_num,N);
%             for ii = 1:N+1
%                 for jj = 1:this.gmm_num
%                     P_pred{jj,ii} = sdpvar(2,2,'full'); % a set of 2-by-2 symmetric matrices
%                 end
%             end
            
            % auxiliary variable
            K = sdpvar(2*this.gmm_num,2*N,'full');
            
%             while (1)
                
                % obj
                obj = 0;%P(1,1,N+1)+P(2,2,N+1); % trace of last covariance
                for ii = 1:N
                    for jj = 1:this.gmm_num
                        tmp = 0;
                        for ll = 1:this.gmm_num
                            % obj uses the 0-order approximation
                            %%% I assume that the covariance does not change
                            %%% for now, which is the simplification. Will
                            %%% change this later after making program work.
                            if ll == jj
                                tmp = tmp+this.wt(ll)/(2*pi*det(this.P{ll}));
                            else
                                tmp = tmp+this.wt(ll)/(2*pi*det(this.P{ll}))*exp(-(x(2*jj-1:2*jj,ii+1)...
                                    -x(2*ll-1:2*ll,ii+1))'/this.P{ll}*(x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))/2);
                            end
                        end
                        obj = obj-this.wt(jj)*log(tmp); % missing term: 1/2*E((x-mu)^T*g''(mu)*(x-mu))
                    end
                end
                
                % constraints
                % initial value
                constr = [[z(:,1) == this.state]:'init_z'];
                constr = [constr,[x(:,1) == this.est_pos(:)]:'init_x'];
                for jj = 1:this.gmm_num
                    constr = [constr,[P(:,:,jj,1) == this.P{jj}]:'init_P'];%[1 0;0 1]];
%                     constr = [constr,[P(2*jj-1:2*jj,1:2) == this.P{jj}]:'init_P'];%[1 0;0 1]];
%                     constr = [constr,[P{jj,1} == this.P{jj}]:'init_P'];%[1 0;0 1]];
                end
                
                % constraints on the go
                for ii = 1:N
                    % robot state
%                     if isempty(zref)
                        constr = [constr,[z(:,ii+1) == z(:,ii)+...
                            [z(4,ii)*cos(z(3,ii));z(4,ii)*sin(z(3,ii));...
                            u(:,ii)]*dt]:'robot motion'];
%                     else
%                         % linearize using previous result
%                         constr = [constr,z(:,ii+1) == z(:,ii)+...
%                             [z(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(z(3,ii)-zref(3,ii));
%                             z(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(z(3,ii)-zref(3,ii));
%                             u(:,ii)]*dt];
%                     end
                    
                    constr = [constr,[[fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=...
                        [fld.fld_cor(2);fld.fld_cor(4)]]:'state bound'];                                       
                    
                    % target prediction
                    for jj = 1:this.gmm_num
                        A = del_f(x(2*jj-1:2*jj,ii));
%                         if isempty (zref)
                            C = del_h(x(2*jj-1:2*jj,ii),z(1:2,ii));
%                         else
%                             C = del_h(x(2*jj-1:2*jj,ii),zref(1:2,ii));
%                         end
                        
                        % forward prediction
                        % mean
                        constr = [constr, [x_pred(2*jj-1:2*jj,ii) == f(x(2*jj-1:2*jj,ii))]:'pred_mean'];
                        %                         x_pred = f(x(2*jj-1:2*jj,ii));
                        % covariance
                        constr = [constr,[P_pred(:,:,jj,ii) == A*(P(:,:,jj,ii))*A'+Q]:'pred_cov'];
                        %                         constr = [constr,[P_pred(2*jj-1:2*jj,2*ii-1:2*ii) == A*(P(2*jj-1:2*jj,2*ii-1:2*ii))*A'+Q]:'pred_cov'];
%                         constr = [constr,[P_pred{jj,ii} == A*P{jj,ii}*A'+Q]:'pred_cov'];
                        %                         P_pred = A*P{jj,ii}*A'+Q;
                        
                        % update using pesudo measurement
                        %                         T = C*P_pred*C'+R;
                        T = C*P_pred(:,:,jj,ii)*C'+R;
%                         T = C*P_pred{jj,ii}*C'+R;
                        constr = [constr, [K(2*jj-1:2*jj,2*ii-1:2*ii)*T == P_pred(:,:,jj,ii)*C']:'K'];
                        %                         constr = [constr, [K(2*jj-1:2*jj,2*ii-1:2*ii)*T == P_pred(2*jj-1:2*jj,2*ii-1:2*ii)*C']:'K']; % define K=P_pred*C'(C*P_pred*C'+T)^-1
%                         constr = [constr, [K(2*jj-1:2*jj,2*ii-1:2*ii)*T == P_pred{jj,ii}*C']:'K'];
                        
                        
                        % mean
                        %                         constr = [constr,[x(2*jj-1:2*jj,ii+1) == x_pred(2*jj-1:2*jj,ii)]:'upd_mean'];
                        %                         constr = [constr,[x(2*jj-1:2*jj,ii+1) == x_pred]:'upd_mean'];
                        constr = [constr,[x(2*jj-1:2*jj,ii+1) == x_pred(2*jj-1:2*jj,ii)]:'upd_mean'];                        
                        
                        % covariance
                        % if the robot can replan reasonably fast, then it
                        % can be fine to use the initial theta_bar as the
                        % approximation
                        
%                         if isempty(zref) 
                            constr = [constr,[(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii))*gam_den(z(1:2,ii+1),z(3,ii+1),...
                                x_olp(2*jj-1:2*jj,ii+1),this.alp1,this.alp2,this.alp3)...
                                == -K(2*jj-1:2*jj,2*ii-1:2*ii)*C*P_pred(:,:,jj,ii)]:'upd_cov'];
%                             constr = [constr,[(P{jj,ii+1}-P_pred{jj,ii})*gam_den(z(1:2,ii+1),z(3,ii+1),...
%                                 x_ol_pred(2*jj-1:2*jj,ii+1),theta_bar(jj),this.alp1,this.alp2)...
%                                 == -K(2*jj-1:2*jj,2*ii-1:2*ii)*C*P_pred{jj,ii}]:'upd_cov'];
%                         end
                    end
                end
                constr = [constr, [this.w_lb <= u(1,:) <= this.w_ub, this.a_lb <= u(2,:) <= this.a_ub...
                    this.v_lb <= z(4,:) <= this.v_ub]:'input_bound'];
                
                
                % seed solution for current iteration
                assign(z,init_sol.z);
                assign(u,init_sol.u);
                assign(x,init_sol.x);
                assign(K,init_sol.K);
                assign(P,init_sol.P);
                assign(P_pred,init_sol.P_pred);
                assign(x_pred,init_sol.x_pred);
                
%                 if ~isempty(zref)
%                     assign(z,zref)
%                     assign(u,uref)
%                 end
                opt = sdpsettings('solver','ipopt','verbose',3,'debug',1,'showprogress',1,'usex0',1);
                
                sol1 = optimize(constr,obj,opt);
                zref = value(z);
                uref = value(u);
                Kref = value(K);
                xref = value(x);
                P_pred = value(P_pred);

%                 dif = norm(is_in_fov-is_in_fov_approx,1);
%                 alp = alp*alp_inc;
                
                % check the singularity of P
%                 for ii = 1:N
%                     for jj = 1:this.gmm_num
%                         display(cond(value(P{jj,ii})))
%                     end
%                 end
                
                display('Robot.m line 780')
                display(value(P))
                
                init_sol = struct('zref',zref,'uref',uref,'Kref',Kref,'xref',xref,'P_pred_ref',P_pred);
                [optz,optu] = cvxPlanner(this,fld,init_sol);
                
%             end
        end
        
        function [init_state] = genInitState(this,fld,prev_state)
            % generate initial states for cvxPlanner
            % quantities to compute: x,z,u,K,P,P_pred,x_pred, theta_bar
            N = this.mpc_hor;
            
            tar = fld.target;
            Q = tar.Q;
            del_f = tar.del_f;
            f = tar.f;
            R = this.R;
            del_h = this.del_h;
            dt = this.dt;
            
            % open-loop prediction of target position. no measurement info 
            % is used
            % get x, x_pred
            x_olp = zeros(2*this.gmm_num,N+1);
            x_olp(:,1) = this.est_pos(:);
            for ii = 1:N
                x_olp(:,ii+1) = f(x_olp(:,ii));
            end            
            init_state.x_olp = x_olp;
            init_state.x = x_olp;
            init_state.x_pred = x_olp(:,2:end);                        
            
            % get z,u
            if isempty(prev_state)
                % if this is the 1st time of calling ngPlanner, use motion
                % primitives
                
                u = [0;0.5]*ones(1,N);
                z = zeros(4,N+1);
                
                z(:,1) = this.state;
                
                for ii = 1:N
                    z(:,ii+1) = z(:,ii)+...
                        [z(4,ii)*cos(z(3,ii));z(4,ii)*sin(z(3,ii));...
                        u(:,ii)]*dt;
                end
            else
                
                optz = prev_state.optz;
                optu = prev_state.optu;           
                % if there exist previous solution, extrapolate them
                u = [optu(:,2:end),optu(:,end)];
                tmp_z = optz(:,end)+[optz(4,end)*cos(optz(3,end));optz(4,end)*sin(optz(3,end));...
                    optu(:,end)]*dt;
                z = [optz(:,2:end),tmp_z];
            end
            init_state.z = z;
            init_state.u = u;
            
            % get theta_bar
            theta_bar = zeros(this.gmm_num,N+1);
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    tmp_vec = x_olp(2*jj-1:2*jj,ii)-z(1:2,ii);
                    theta_bar(jj,ii) = atan2(tmp_vec(2),tmp_vec(1));
                end
            end
            init_state.theta_bar = theta_bar;
            
            % get P, P_pred, K
            P = zeros(2,2,this.gmm_num,N+1);
            P_pred = zeros(2,2,this.gmm_num,N);
            K = zeros(2*this.gmm_num,2*N);
            
            for jj = 1:this.gmm_num
                P(:,:,jj,1) = this.P{jj};
            end
            
            for ii = 1:N                
                for jj = 1:this.gmm_num
                    A = del_f(x_olp(2*jj-1:2*jj,ii+1));
                    C = del_h(x_olp(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    P_pred(:,:,jj,ii) = A*P(:,:,jj,ii)*A'+Q;
                    T = C*P_pred(:,:,jj,ii)*C'+R;
                    K(2*jj-1:2*jj,2*ii-1:2*ii)= P_pred(:,:,jj,ii)*C'/T;
                    gam = this.inFOV(x_olp(2*jj-1:2*jj,ii));
                    P(:,:,jj,ii+1) = P_pred(:,:,jj,ii)-gam*K(2*jj-1:2*jj,2*ii-1:2*ii)*P_pred(:,:,jj,ii)*C;
                end
            end
            init_state.P = P;
            init_state.P_pred = P_pred;
            init_state.K = K;       
            
            % visualize the seed path
            hold on
            this.plotSeedTraj(z,x_olp,fld)
        end
               
        function [init_state] = genInitState_kf(this,fld,prev_state)
            % generate initial states for cvxPlanner2
            % quantities to compute: x,z,u,K,P,P_pred,x_pred, theta_bar
            N = this.mpc_hor;
            
            tar = fld.target;
            Q = tar.Q;
            del_f = tar.del_f;
            f = tar.f;
            R = this.R;
            del_h = this.del_h;
            dt = this.dt;
            
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            
            % open-loop prediction of target position. no measurement info 
            % is used
            % get x, x_pred
            x_olp = zeros(2,N+1);
            x_olp(:,1) = this.est_pos;
            for ii = 1:N
                x_olp(:,ii+1) = f(x_olp(:,ii));
            end            
            init_state.x_olp = x_olp;
            init_state.x = x_olp;
            init_state.x_pred = x_olp(:,2:end);                        
            
            % get z,u
            if isempty(prev_state)
                % if this is the 1st time of calling ngPlanner, use motion
                % primitives
                
                u = [0;1]*ones(1,N);
                z = zeros(4,N+1);
                
                z(:,1) = this.state;
                
                for ii = 1:N
                    z(:,ii+1) = z(:,ii)+...
                        [z(4,ii)*cos(z(3,ii));z(4,ii)*sin(z(3,ii));...
                        u(:,ii)]*dt;
                end
            else
                optz = prev_state.optz;
                optu = prev_state.optu;           
                % if there exist previous solution, extrapolate them
                u = [optu(:,2:end),optu(:,end)];
                tmp_z = optz(:,end)+[optz(4,end)*cos(optz(3,end));optz(4,end)*sin(optz(3,end));...
                    optu(:,end)]*dt;
                z = [optz(:,2:end),tmp_z];
            end
            init_state.z = z;
            init_state.u = u;            
            
            % get P, P_pred, K
            P = zeros(2,2,N+1);
            P_pred = zeros(2,2,N);
            K = zeros(2,2*N);
            
            for jj = 1:this.gmm_num
                P(:,:,1) = this.P;
            end
            
            for ii = 1:N
                A = del_f(x_olp(:,ii+1));
                C = del_h(x_olp(:,ii+1),z(1:2,ii+1));
                P_pred(:,:,ii) = A*P(:,:,ii)*A'+Q;
                T = C*P_pred(:,:,ii)*C'+R;
                K(:,2*ii-1:2*ii)= P_pred(:,:,ii)*C'/T;
%                 gam = this.inFOV(x_olp(:,ii));
                gam = this.gam(z(1:2,ii+1),z(3,ii+1),x_olp(:,ii+1),alp1,alp2,alp3);
                P(:,:,ii+1) = P_pred(:,:,ii)-gam*K(:,2*ii-1:2*ii)*P_pred(:,:,ii)*C;
            end
            init_state.P = P;
            init_state.P_pred = P_pred;
            init_state.K = K;            
            
            % for debug purpose only
            % compute the gamma and gamma_grad for the predicted state
            %{
            fov_ref = zeros(N+1,1);
            gam_ref = zeros(N+1,1);
            gam_grad_ref = zeros(3,N+1);
            for ii = 1:N+1
                tmp = this;
                tmp.state = z(1:3,ii);
                fov_ref(ii) = inFOV(tmp,x_olp(:,ii));
                gam_ref(ii) = gam(this,z(1:2,ii),z(3,ii),x_olp(:,ii),alp1,alp2,alp3);                
                gam_grad_ref(:,ii) = gam_grad(this,z(1:2,ii),z(3,ii),x_olp(:,ii),alp1,alp2,alp3);
            end
            %}
            
            % visualize the seed path
            hold on
            this.plotSeedTraj(z,x_olp,fld)
        end
        
        %% robot state updating
        function this = updState(this,u)
            st = this.state;
            this.optu = [this.optu,u(:,1)];
            dt = this.dt;
            this.state = st+[st(4)*cos(st(3));st(4)*sin(st(3));u(:,1)]*dt;
            this.traj = [this.traj,this.state];
        end               
        
        function z = simState(this,u,z0)
            if nargin > 2
                st = z0;
            else
                st = this.state;
            end
            
            dt = this.dt;
            len = size(u,2);
            z = zeros(length(st),len+1);
            z(:,1) = st;
            for ii = 1:len
                z(:,ii+1) = z(:,ii)+[z(4,ii)*cos(z(3,ii));z(4,ii)*sin(z(3,ii));u(:,ii)]*dt;
            end
        end
        
        %% %%%%% utilities for numerical optimization
        %% scp solver (adapted from CS 287 class)
        function [x, trust_box_size, success, merit, model_merit, new_merit] = minimize_merit_function(this, x, Q, q, f, A_ineq, b_ineq, A_eq, b_eq, glin, hlin, g, h, hinge, abssum, cfg, penalty_coeff, trust_box_size, snum, fld)
            % f: nonlinear, non-quad obj, g: nonlinear inequality constr, h: nonlinear equality constr.
            % glin: linear inequality constr, hlin: linear equality constr.
            
            dim_x = length(x);
            success = true;
            sqp_iter = 1;
            N = this.mpc_hor;
%             del_f = cfg.del_f;
%             del_h = cfg.del_h;
            
%             fquadlin = @(x) q*x + .5*x'*(Q*x);
            %%%%% remember to uncomment the following line when doing NGP
            fquadlin = @(x) this.getQuadObj(x,snum); 
%             hinge = @(x) sum(max(x,0));
%             abssum = @(x) sum(abs(x));
            
            while  true
                % In this loop, we repeatedly construct a quadratic approximation
                % to the nonlinear part of the objective f and a linear approximation to the nonlinear
                % constraints f and g.
                fprintf('      [sqp loop] Robot.m, line %d\n', MFileLineNr())
                fprintf('      sqp iter: %i\n', sqp_iter);

                if cfg.f_use_numerical % && ~isempty(f)
                    fval = f(x);
                    [fgrad, fhess] = this.getNumGradHess(f,x,cfg.full_hessian);
                    % diagonal adjustment
                    mineig = min(eigs(fhess));
%                     if mineig < 0
%                         fprintf('[sqp loop] Robot.m, line %d', MFileLineNr())
%                         fprintf('negative hessian detected. adjusting by %.3g\n',-mineig);
%                         fhess = fhess + eye(dim_x) * ( - mineig);
%                     end
                else %if ~isempty(f)
                    [fval, fgrad, fhess] = f(x);
%                 else
%                     fval = 0;
%                     fgrad = 0; 
%                     fhess = 0;
                end
                
                if cfg.g_use_numerical% && ~isempty(g)
                    gval = g(x);
                    gjac = this.getNumJac(g,x);
                else %if ~isempty(g)
                    [gval, gjac] = g(x);
%                 else
%                     gval = 0;
%                     gjac = 0;
                end
                
                if cfg.h_use_numerical% && ~isempty(h)
                    hval = h(x);
                    hjac = this.getNumJac(h,x);
                else %if ~isempty(h)
                    [hval, hjac] = h(x);
%                 else
%                     hval = 0;
%                     hjac = 0;
                end
                
                merit = fval + fquadlin(x) + penalty_coeff * ( hinge(gval) + abssum(hval) );
                
                tmp_cnt = 1; % a temp counter for debugging use only
                
                while true
                    % This is the trust region loop
                    % Using the approximations computed above, this loop shrinks
                    % the trust region until the progress on the approximate merit
                    % function is a sufficiently large fraction of the progress on
                    % the exact merit function.
                    
                    fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('        trust region size: %.3g\n', trust_box_size);
                    fprintf('        iteration %d with current sqp loop\n', tmp_cnt);
                    
                    % YOUR CODE INSIDE CVX_BEGIN and CVX_END BELOW
                    % Write CVX code to minimize the convex approximation to
                    % the merit function, using the jacobians computed above.
                    % It should create variable xp, which is the candidate for
                    % updating x -> xp.
                    
                    % You should enforce the linear constraints exactly.
                    % Make sure to include the constant term f(x) in the merit function
                    % objective as the resulting cvx_optval is used further below.
                    
                    
                    cvx_begin quiet %YOUR CODE HERE
                        variables xp(dim_x, 1);
                        %                 minimize fquadlin(x) + fval + (x'*Q'+q+fgrad)*(xp-x) + penalty_coeff*( hinge(gval+gjac*(xp-x)) + abssum(hval+hjac*(xp-x)) )
                        minimize fquadlin(xp) + fval + (fgrad)*(xp-x) + penalty_coeff*( hinge(gval+gjac*(xp-x)) + abssum(hval+hjac*(xp-x)) )
                        subject to

                        %                     A_ineq*xp <= b_ineq;
                        %                     A_eq*xp == b_eq;
                        hlin(xp) == 0;
                        glin(xp) <= 0;
                        norm(xp-x,2) <= trust_box_size;
                    cvx_end
                    
                    tmp_cnt = tmp_cnt+1;
                    
                    if strcmp(cvx_status,'Failed')
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        Failed to solve QP subproblem.\n');
                        success = false;
                        return;
                    elseif contains(cvx_status,'Infeasible')
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        Infeasible QP subproblem occurred.\n');
                        
                    end
                    
                    fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('        cvx status: %s\n',cvx_status)
                    
                    model_merit = cvx_optval;
                    new_merit = f(xp) + fquadlin(xp) + penalty_coeff * ( hinge(g(xp)) + abssum(h(xp)) ) ;
                    approx_merit_improve = merit - model_merit;
                    exact_merit_improve = merit - new_merit;
                    merit_improve_ratio = exact_merit_improve / approx_merit_improve;
                    
%                     sprintf('[inner sqp loop] Robot.m, line %d', MFileLineNr())
%                     info = struct('trust_box_size',trust_box_size);
                    
                    fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('        approx improve: %.3g. exact improve: %.3g. ratio: %.3g\n', approx_merit_improve, exact_merit_improve, merit_improve_ratio);
                    if approx_merit_improve < -1e-5
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        Approximate merit function got worse (%.3e).\n',approx_merit_improve);
                        fprintf('        Either convexification is wrong to zeroth order, or you''re in numerical trouble\n');
                        fprintf('        Return to penalty loop\n');
                        success = false;
                        return;
                    elseif approx_merit_improve < cfg.min_approx_improve
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        Converged: y tolerance\n');
                        fprintf('        Return to penalty loop\n');
                          
                        x = xp;
%                         if ~isempty(cfg.callback), cfg.callback(x,info); end
                        return;
                    elseif (exact_merit_improve < 0) || (merit_improve_ratio < cfg.improve_ratio_threshold)
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        shrink trust region. stay in current sqp iteration\n');
                        trust_box_size = trust_box_size * cfg.trust_shrink_ratio;
                    else
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        expand trust region. go to next sqp iteration\n');
                        if abs(norm(x-xp,2) - trust_box_size) < 1e-3
                            trust_box_size = trust_box_size * cfg.trust_expand_ratio;
                        end
                        x = xp;
%                         if ~isempty(cfg.callback), cfg.callback(x,info); end
                        break; % from trust region loop
                    end
                    
                    
                    if trust_box_size < cfg.min_trust_box_size
                        fprintf('        [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
                        fprintf('        Converged: x tolerance\n');
                        fprintf('        Return to penalty loop\n');
                        return;
                    end
                    
                    % visualize solution
                    
%                     this.plotPlannedTraj(this.convState(x,snum,'z'),[],fld)
                    
                end % tr
                
                sqp_iter = sqp_iter + 1;
                
                % update state by using computed u
%                 x2 = this.updOptVar(x,snum);
%                 sprintf('after updating, the different between x2 and x is %d',norm(x2-x,2))
%                 x = x2;
               
                if sqp_iter > cfg.max_sqp_iter
                    break
                end
            end % sqp
        end
           
        function [x,success] = find_closest_feasible_point(this, x0, glin, hlin)
            % Find a point that satisfies linear constraints, if x0 doesn't
            % satisfy them
            
            success = true;
            if any(glin(x0) > 0) || any(hlin(x0) ~= 0)
                fprintf('Robot.m, line %d\n', MFileLineNr())
                fprintf('initialization doesn''t satisfy linear constraints. finding the closest feasible point\n');
                
                cvx_begin quiet
                variables('x(length(x0))')
                minimize('sum((x-x0).^2)')
                subject to
                    hlin(x) == 0;
                    glin(x) <= 0;
                cvx_end
                
                if strcmp(cvx_status,'Failed') || strcmp(cvx_status,'Infeasible')
                    success = false;
                    fprintf('Robot.m, line %d\n', MFileLineNr())
                    fprintf('Couldn''t find a point satisfying linear constraints\n')
                    return;
                end
            else
                x = x0;
            end
        end
        
        % update optimization variables for another loop using open-loop
        % updating
        function snew = updOptVar(this,s,snum,alp1,alp2,alp3)  
            
            z_orig = this.convState(s,snum,'z');
            u_orig = this.convState(s,snum,'u');
            x_orig = this.convState(s,snum,'x'); %s(snum(3,1),snum(3,2));
%             xpred = this.convState(s,snum,'x_pred'); %s(snum(5,1),snum(5,2));
            P_orig = this.convState(s,snum,'P'); %s(snum(3,1),snum(3,2));
%             P_pred = this.convState(s,snum,'P_pred'); %s(snum(5,1),snum(5,2));
            
            N = this.mpc_hor;
            
%             alp1 = this.alp1;
%             alp2 = this.alp2;
%             alp3 = this.alp3;
            del_h = this.del_h;
            R = this.R;
            f = this.target.f;
            del_f = this.target.del_f;
            Q = this.target.Q;

            % z
            z = this.simState(u_orig,z_orig(:,1));
            
            % belief state
            x = zeros(2,N+1);
            P = zeros(2,2,N+1);
            P_pred = zeros(2,2,N);
            K = zeros(2,2*N);
            
            % initialization
            x(:,1) = x_orig(:,1);
            P(:,:,1) = P_orig(:,:,1);
            
            for ii = 1:N
                x(:,ii+1) = f(x(:,ii));
            end
            
            xpred = x(:,2:end);
            
            for ii = 1:N
                A = del_f(x(:,ii+1));
                C = del_h(x(:,ii+1),z(1:2,ii+1));
                P_pred(:,:,ii) = A*P(:,:,ii)*A'+Q;
                T = C*P_pred(:,:,ii)*C'+R;
                K(:,2*ii-1:2*ii)= P_pred(:,:,ii)*C'/T;
%                 gam = this.inFOV(x_olp(:,ii));
                gam = this.gam(z(1:2,ii+1),z(3,ii+1),x(:,ii+1),alp1,alp2,alp3);
                P(:,:,ii+1) = P_pred(:,:,ii)-gam*K(:,2*ii-1:2*ii)*C*P_pred(:,:,ii);
            end
            snew = setState(this,z,u_orig,x,xpred,P,P_pred,K);
        end
        
        %% objective function
        % for general NGP. this one may duplicate cmpObj. the difference is 
        % mainly on the parameters.
        function val = getObj(this,s,snum,alp1,alp2,alp3)
            x = this.convState(s,snum,'x'); %s(snum(3,1),snum(3,2));
            u = this.convState(s,snum,'u'); %s(snum(3,1),snum(3,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1),snum(5,2));
            z = this.convState(s,snum,'z');
            
            N = this.mpc_hor;
            val = 0;
            for ii = 2:N+1
               for jj = 1:this.gmm_num
                   tmp = 0;
                   for ll = 1:this.gmm_num                       
                       P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/2; % numerical issues happen that makes P non symmetric
%                        sprintf('Robot.m, line %d', MFileLineNr())
%                        display(P(:,:,ll,ii))
                       mineigv = min(eig(P(:,:,ll,ii)));                       
                       if mineigv <= 0
                           P(:,:,ll,ii) = P(:,:,ll,ii)+(-mineigv+0.01)*eye(size(P(:,:,ll,ii),1));
                       end                       
                       tmp = tmp+mvnpdf(x(2*jj-1:2*jj,ii),x(2*ll-1:2*ll,ii),P(:,:,ll,ii));
                   end
%                    tmp_dis = sum((x(:,ii)-z(1:2,ii)).^2);
                   tmp_dis = abs(sum((x(:,ii)-z(1:2,ii)).^2)-1); % distance between sensor and MLE target position
                   
                   % added on 8/20
                   % temporarily add the diff between gamma and 1 as
                   % penalty
                   tar_pos = x(2*jj-1:2*jj,ii);                  

                   val = val-this.wt(jj)*log(tmp)+tmp_dis+(1-this.gam(z(1:2,ii),z(3,ii),tar_pos,alp1,alp2,alp3));
               end               
               val = val+sum(u(:,ii-1).^2); % penalize on control input
            end
            val = 0.1*val;
%             val = 0;
        end    
        
        % for KF
        % not in use for KF case now
        function val = getObj_kf(this,s,snum) % not in use
            x = this.convState(s,snum,'x'); %s(snum(3,1),snum(3,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1),snum(5,2));
            
            N = this.mpc_hor;
            val = 0;
            for ii = 1:N+1
               for jj = 1:this.gmm_num
                   val = val+trace(P(:,:,jj,ii));
               end
            end
            val = -val;
        end    
        
        function val = getQuadObj(this,x,snum)
            % linear quadratic term in the objective function
            val = 0;
%             u = this.convState(x,snum,'u');
%             val = sum(sum(u.^2));
        end
        
        % compute the exact value (using 0-th approx) of objective function
        function val = cmpObj(this,z,P)
            N = this.mpc_hor;
            
            if (nargin < 3)
                % this case means x is an augmented parameter containing x and P              
                x = reshape(z(1:2*this.gmm_num*(N+1)),2*this.gmm_num,N+1);
                P = reshape(z(2*this.gmm_num*(N+1)+1:end),2,2,this.gmm_num,N+1);                
            else 
                x = z;
            end
            
            val = 0;
            for ii = 2:N+1
               for jj = 1:this.gmm_num
                   tmp = 0;
                   for ll = 1:this.gmm_num
                       P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/2; % numerical issues happen that makes P non symmetric
                       tmp = tmp+mvnpdf(x(2*jj-1:2*jj,ii),x(2*ll-1:2*ll,ii),P(:,:,ll,ii));
                   end                   
                   val = val+this.wt(jj)*log(tmp);
               end
            end
            val = -val;
        end
        
        % compute the approximate value of objective function using the
        % same formula as the cvx obj, which is the upper found of cmpObj
        function val = cmpObjAprx(this,x,P)
            N = this.mpc_hor;
            val = 0;
            for ii = 2:N+1
               for jj = 1:this.gmm_num    
                   tmp = 0;
                   for ll = 1:this.gmm_num
%                        tmp = tmp+log(mvnpdf(x(2*ll-1:2*ll,ii),x(2*jj-1:2*jj,ii),P(:,:,jj,ii)));
                       tmp = tmp+log(mvnpdf(x(2*jj-1:2*jj,ii),x(2*ll-1:2*ll,ii),P(:,:,ll,ii)));
                       display(P(:,:,ll,ii))
                   end
                   val = val+this.wt(jj)*tmp;
               end
            end        
            val = -val;
        end
        
        % compute the approximate merit function which put nonlinear
        % equality and inequalities into the objective function (but not
        % linearizing them)
         function val = cmpMerit(this,z,u,zref,x,P,P_pred,Kref,mu)
            N = this.mpc_hor;
            val = 0;
            for ii = 2:N+1
               for jj = 1:this.gmm_num    
                   tmp = 0;
                   for ll = 1:this.gmm_num
                       P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/ 2; % numerical issues happen that makes P non symmetric
                       tmp = tmp+mvnpdf(x(2*jj-1:2*jj,ii),x(2*ll-1:2*ll,ii),P(:,:,ll,ii));
                   end
                   val = val+this.wt(jj)*log(tmp);
               end
            end        
            val = -val;
            val = val + mu*(this.penKinConstr(z,u,zref)+this.penBelConstr(z,x,P,P_pred,Kref));
         end
        
        %% constraints
        % for general NGP
        function h = getKinConstr(this,s,snum)
            % kinematic constraint (nonlinear equality)
            z = this.convState(s,snum,'z'); %s(snum(1,1):snum(1,2));
            u = this.convState(s,snum,'u'); %s(snum(2,1):snum(2,2));
            h = [];
            N = this.mpc_hor;
            dt = this.dt; 
            for ii = 1:N
%                 val = val+z(:,ii+1) - (z(:,ii)+ [z(4,ii)*cos(z(3,ii)) 0; z(4,ii)*sin(z(3,ii)) 0;...
%                     1 0; 0 1]*u(:,ii)*dt);
                h = [h;z(:,ii+1) - (z(:,ii)+ [z(4,ii)*cos(z(3,ii)); z(4,ii)*sin(z(3,ii));...
                    u(:,ii)]*dt)];
            end
        end
        
        function h = getBelConstr(this,s,snum,Kref,alp1,alp2,alp3)
            % belief update constraint (nonlinear equality)
            z = this.convState(s,snum,'z'); %s(snum(1,1):snum(1,2));
            u = this.convState(s,snum,'u'); %s(snum(2,1):snum(2,2));
            x = this.convState(s,snum,'x'); %s(snum(3,1):snum(3,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1):snum(5,2));
            P_pred = this.convState(s,snum,'P_pred'); %s(snum(6,1):snum(6,2));
            h = 0;
            
            N = this.mpc_hor;            
            dt = this.dt;
            gam = @this.gam;
            del_h = this.del_h;
%             alp1 = this.alp1;
%             alp2 = this.alp2;
%             alp3 = this.alp3;
            
            % xpred_k+1 = f(x_k)
            %%% this part is now implemened as a linear equality constraint
            %%% will be added later when considering a moving target
            
            % P_k+1|k+1=P_k+1|k-gamma*K*C*P_k+1|k
            C = zeros(2,2,this.gmm_num,N);
            for ii = 1:N
                % target prediction
                for jj = 1:this.gmm_num
                    C(:,:,jj,ii) = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    
                    % covariance
                    tmp_sum = zeros(2,2);
                    T = Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C(:,:,jj,ii);
                    
                    for ll = 1:this.gmm_num
                        %%% note: gamma depends on ll, C depends
                        %%% only on jj
                        tar_pos = x(2*ll-1:2*ll,ii+1); % use each gmm component mean as a possible target position                        
                        tmp_sum = tmp_sum+this.wt(ll)*gam(z(1:2,ii+1),z(3,ii+1),...
                        tar_pos,alp1,alp2,alp3)*T*P_pred(:,:,jj,ii);
                    end
                    h = h+sum(sum(abs(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum)));
                end
            end
        end
        
        function hlin = getLinEqConstr(this,s,snum,tar)
            % linear equality constraints
            z = this.convState(s,snum,'z'); %s(snum(1,1):snum(1,2));
            x = this.convState(s,snum,'x'); %s(snum(3,1):snum(3,2));
            x_pred = this.convState(s,snum,'x_pred'); %s(snum(4,1):snum(4,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1):snum(5,2));
            P_pred = this.convState(s,snum,'P_pred'); %s(snum(6,1):snum(6,2));            
            
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            N = this.mpc_hor;
            hlin = [];   
            % initial condition
            % z(:,1) == this.state;
            % x(:,1) == this.est_pos(:);
            % for jj = 1:this.gmm_num
              %  triu(P(:,:,jj,1)) == triu(this.P{jj});
            % end
            hlin = [hlin;z(:,1)-this.state];
            hlin = [hlin;x(:,1)-this.est_pos(:)];
            for jj = 1:this.gmm_num
               tmp = triu(P(:,:,jj,1))-triu(this.P{jj});
               hlin = [hlin;tmp(:)];
            end
            
            A = zeros(2,2,this.gmm_num,N);
            for ii = 1:N
                for jj = 1:this.gmm_num                    
                    % forward prediction
                    % x_k+1|k = f(x_k)
                    hlin = [hlin;x_pred(2*jj-1:2*jj,ii)-f(x(2*jj-1:2*jj,ii))];
                    % P_k+1|k=A*P_k*A+Q
                    A(:,:,jj,ii) = del_f(x(2*jj-1:2*jj,ii));
                    tmp = triu(P_pred(:,:,jj,ii))-triu(A(:,:,jj,ii)*P(:,:,jj,ii)*A(:,:,jj,ii)'+Q);%+triu(slk_P_pred(:,:,jj,ii));
                    hlin = [hlin;tmp(:)];
                    % update
                    % x_k+1|k+1 = x_k+1|k
                    hlin = [hlin;x(2*jj-1:2*jj,ii+1)-x_pred(2*jj-1:2*jj,ii)];
                end
            end
        end
        
        function glin = getLinIneqConstr(this,s,snum,fld)
            % linear inequality constraints
            z = this.convState(s,snum,'z'); %s(snum(1,1),snum(1,2));
            u = this.convState(s,snum,'u');%s(snum(2,1),snum(2,2));
            
            N = this.mpc_hor;            
%             glin = [];
            % bounds on states and input
            % this.w_lb <= u(1,:) <= this.w_ub;
            % this.a_lb <= u(2,:) <= this.a_ub;
            % this.v_lb <= z(4,:) <= this.v_ub;            
            glin = [(u(1,:)-this.w_ub)';(u(2,:)-this.a_ub)';(z(4,:)-this.v_ub)';...
                (this.w_lb-u(1,:))';(this.a_lb-u(2,:))';(this.v_lb-z(4,:))'];
            % [fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=[fld.fld_cor(2);fld.fld_cor(4)];
            for ii = 1:N
                glin = [glin;z(1:2,ii+1)-[fld.fld_cor(2);fld.fld_cor(4)];...
                    [fld.fld_cor(1);fld.fld_cor(3)]-z(1:2,ii+1)];
            end
        end
        
        % for KF        
        function val = getLinEqConstr_kf(this,s,snum,tar)
            % linear equality constraints for KF
            z = this.convState(s,snum,'z'); %s(snum(1,1):snum(1,2));
            x = this.convState(s,snum,'x'); %s(snum(3,1):snum(3,2));
            x_pred = this.convState(s,snum,'x_pred'); %s(snum(4,1):snum(4,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1):snum(5,2));
            P_pred = this.convState(s,snum,'P_pred'); %s(snum(6,1):snum(6,2));            
            
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            N = this.mpc_hor;
            val = [];
            val1 = [];
            val2 = [];
            val3 = [];
            % initial condition
            % z(:,1) == this.state;
            % x(:,1) == this.est_pos(:);
            % for jj = 1:this.gmm_num
              %  triu(P(:,:,jj,1)) == triu(this.P{jj});
            % end
            val = [val;z(:,1)-this.state];
            val = [val;x(:,1)-this.est_pos(:)];
            for jj = 1:this.gmm_num
               tmp = triu(P(:,:,jj,1))-triu(this.P); %%% note: triu here returns a 2*2 matrix, with lower-left item being 0
               val = [val;tmp(:)];
            end
            
            A = zeros(2,2,this.gmm_num,N);
            for ii = 1:N
                for jj = 1:this.gmm_num                    
                    % forward prediction
                    % x_k+1|k = f(x_k)
                    val1 = [val1;x_pred(2*jj-1:2*jj,ii)-f(x(2*jj-1:2*jj,ii))];
                    % P_k+1|k=A*P_k*A+Q
                    A(:,:,jj,ii) = del_f(x(2*jj-1:2*jj,ii));
                    tmp = triu(P_pred(:,:,jj,ii))-triu(A(:,:,jj,ii)*P(:,:,jj,ii)*A(:,:,jj,ii)'+Q);%+triu(slk_P_pred(:,:,jj,ii));
                    val2 = [val2;tmp(:)];
                    % update
                    % mean
                    % x_k+1|k+1 = x_k+1|k
                    val3 = [val3;x(2*jj-1:2*jj,ii+1)-x_pred(2*jj-1:2*jj,ii)];                   
                end
            end
            val = [val;val1;val2;val3];
        end

        function h = getBelConstr_kf(this,s,snum,alp1,alp2,alp3)
            % belief update constraint (nonlinear equality)
            z = this.convState(s,snum,'z'); %s(snum(1,1):snum(1,2));
            u = this.convState(s,snum,'u'); %s(snum(2,1):snum(2,2));
            x = this.convState(s,snum,'x'); %s(snum(3,1):snum(3,2));
            P = this.convState(s,snum,'P'); %s(snum(5,1):snum(5,2));
            P_pred = this.convState(s,snum,'P_pred'); %s(snum(6,1):snum(6,2));
            K = this.convState(s,snum,'K');
            
            val1 = [];
            val2 = [];
            
            N = this.mpc_hor;            
            dt = this.dt;
            gam = @this.gam;
            del_h = this.del_h;
            R = this.R;
            
            % xpred_k+1 = f(x_k)
            %%% this part is now implemened as a linear equality constraint
            %%% will be added later when considering a moving target
            
            
            C = zeros(2,2,this.gmm_num,N);
            for ii = 1:N
                % target prediction
                for jj = 1:this.gmm_num
                    C(:,:,jj,ii) = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    
                    % update K
                    % K = P_k+1|k*C_k+1'(C_k+1*P_k+1|k*C_k+1'+R)^-1
                    tmp = (K(2*jj-1:2*jj,2*ii-1:2*ii)*(C(:,:,jj,ii)*P_pred(:,:,jj,ii)*C(:,:,jj,ii)'+R)-...
                        P_pred(:,:,jj,ii)*C(:,:,jj,ii)'); % used a scaling factor 
                    val1 = [val1;tmp(:)];
                                        
                    % covariance
                    % P_k+1|k+1=P_k+1|k-gamma*K*C*P_k+1|k
                    tmp_sum = zeros(2,2);
                    T = K(2*jj-1:2*jj,2*ii-1:2*ii)*C(:,:,jj,ii);
                    for ll = 1:this.gmm_num
                        %%% note: gamma depends on ll, C depends
                        %%% only on jj
                        tar_pos = x(2*ll-1:2*ll,ii+1); % use each gmm component mean as a possible target position                        
                        tmp_sum = tmp_sum+this.wt(ll)*gam(z(1:2,ii+1),z(3,ii+1),...
                        tar_pos,alp1,alp2,alp3)*T*P_pred(:,:,jj,ii);
                    end
%                     val2 = val2+sum(sum(abs(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum)));
                    tmp2 = triu(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum);
                    val2 = [val2;tmp2(:)];
                end
            end
            h = [val1;val2];
        end

        % compute the value of nonlinear contraints in sqp obj
        function h = penKinConstr(this,z,u,zref)
            % compute the l1 penalty of contraint
            h = 0;
            N = this.mpc_hor;
            dt = this.dt; 
            for ii = 1:N
                h = h+z(:,ii+1) - (z(:,ii)+ [z(4,ii)*cos(z(3,ii)) 0; z(4,ii)*sin(z(3,ii)) 0;...
                    1 0; 0 1]*u(:,ii)*dt);
            end      
            h = sum(abs(h));
        end
        
        function h = penBelConstr(this,z,x,P,P_pred,Kref)
            % compute the l1 penalty of contraint
            h = 0;
            N = this.mpc_hor;
            dt = this.dt;
            gam = @this.gam;
            del_h = this.del_h;
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            
            for ii = 1:N
                % target prediction
                for jj = 1:this.gmm_num
                    %%%%% note: this part may need change later. In
                    %%%%% fact, linearziation should be wrt
                    %%%%% reference values
                    C = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    
                    % covariance
                    tmp_sum = zeros(2,2);
                    T = Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C;
                    
                    for ll = 1:this.gmm_num
                        %%% note: gamma depends on ll, C depends
                        %%% only on jj
                        tar_pos = x(2*ll-1:2*ll,ii+1); % use each gmm component mean as a possible target position                        
                        tmp_sum = tmp_sum+this.wt(ll)*gam(z(1:2,ii+1),z(3,ii+1),...
                        tar_pos,alp1,alp2,alp3)*T*P_pred(:,:,jj,ii);
                    end
                    
                    h = h+sum(sum(abs(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum)));
                end
            end
        end
        
        %% conversion between the collection of states and seperate states
        function t = convState(this,s,snum,schar)
           % convert the state collection s into a corresponding state variable
           N = this.mpc_hor;
           switch schar
               case 'z'
                   t = s(snum.idx(1,1):snum.idx(1,2));
                   t = reshape(t,4,N+1);                   
               case 'u'
                   t = s(snum.idx(2,1):snum.idx(2,2));
                   t = reshape(t,2,N);
               case 'x'
                   t = s(snum.idx(3,1):snum.idx(3,2));
                   t = reshape(t,2*this.gmm_num,N+1);
               case 'x_pred'
                   t = s(snum.idx(4,1):snum.idx(4,2));
                   t = reshape(t,2*this.gmm_num,N);
               case 'P'
                   t = s(snum.idx(5,1):snum.idx(5,2));
                   t = reshape(t,2,2,this.gmm_num,N+1);
               case 'P_pred'
                   t = s(snum.idx(6,1):snum.idx(6,2));
                   t = reshape(t,2,2,this.gmm_num,N);
               case 'K'
                   t = s(snum.idx(7,1):snum.idx(7,2));
                   t = reshape(t,2,2,this.gmm_num,N);
           end
        end
        
        function [s,snum] = setState(this,z,u,x,xpred,P,P_pred,K)
            % s is the collection of all decision variables
            % s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:)]               
            s = [z(:);u(:);x(:);xpred(:);P(:);P_pred(:);K(:)];
            
            % size of each decision variables
            snum = struct();
            snum.z = length(z(:)); %4*(N+1);
            snum.u = length(u(:)); %2*N;
            snum.x = length(x(:)); %2*this.gmm_num*(N+1);
            snum.xpred = length(xpred(:)); %2*this.gmm_num*N;
            snum.P = length(P(:)); %2*2*this.gmm_num*(N+1);
            snum.P_pred = length(P_pred(:)); %2*2*this.gmm_num*N;
            snum.K = length(K(:)); %2*2*this.gmm_num*N;
            
            % compute the section of s that corresponds to different
            % decision variables
            snum.idx = zeros(7,2);
            snum.idx(1,:) = [1,snum.z];
            snum.idx(2,:) = [snum.idx(1,2)+1,snum.idx(1,2)+snum.u];
            snum.idx(3,:) = [snum.idx(2,2)+1,snum.idx(2,2)+snum.x];
            snum.idx(4,:) = [snum.idx(3,2)+1,snum.idx(3,2)+snum.xpred];
            snum.idx(5,:) = [snum.idx(4,2)+1,snum.idx(4,2)+snum.P];
            snum.idx(6,:) = [snum.idx(5,2)+1,snum.idx(5,2)+snum.P_pred];
            snum.idx(7,:) = [snum.idx(6,2)+1,snum.idx(6,2)+snum.K];
            snum.total = snum.idx(7,2);
        end
        
        %% numerical gradient, Jacobian, and hessian
        function [grad,hess] = getNumGradHess(this,f,s,full_hessian)
            % modified from Pieter Abbeel's CS287
%             z = [x(:);P(:)];
%             f = @this.cmpObj;
            y = f(s);
            assert(length(y)==1);
            
            grad = zeros(1,length(s));
            hess = zeros(length(s));
            
            eps = 1e-5;
            sp = s;
            
            if nargout > 1
                if ~full_hessian
                    for i=1:length(s)
                        sp(i) = s(i) + eps/2;
                        yhi = f(sp);
                        sp(i) = s(i) - eps/2;
                        ylo = f(sp);
                        sp(i) = s(i);
                        hess(i,i) = (yhi + ylo - 2*y)/(eps.^2 / 4);
                        grad(i) = (yhi - ylo) / eps;
                    end
                else
                    grad = this.numerical_jac(f,s);
                    hess = this.numerical_jac(@(z) this.numerical_jac(f,z), s);
                    hess = (hess + hess')/2;
                end
                mineigv = min(eig(hess));
                hess = hess+(mineigv+0.1)*eye(length(s));
            end
        end
        
        % jacobian
        function grad = getNumJac(this,f,s)
            % modified from Pieter Abbeel's CS287            
            y = f(s);
            
            grad = zeros(length(y), length(s));
            
            eps = 1e-5;
            sp = s;
            
            for i=1:length(s)
                sp(i) = s(i) + eps/2;
                yhi = f(sp);
                sp(i) = s(i) - eps/2;
                ylo = f(sp);
                sp(i) = s(i);
                grad(:,i) = (yhi - ylo) / eps;
            end            
        end
        
        % gradient, hessian        
        function [grad,hess] = numerical_grad_hess(this,x,P,full_hessian)
            % modified from Pieter Abbeel's CS287
            z = [x(:);P(:)];
            f = @this.cmpObj;
            y = f(z);
            assert(length(y)==1);
            
            grad = zeros(1, length(z));
            hess = zeros(length(z));
            
            eps = 1e-5;
            zp = z;
            
            if nargout > 1
                if ~full_hessian
                    for i=1:length(z)
                        zp(i) = z(i) + eps/2;
                        yhi = f(zp);
                        zp(i) = z(i) - eps/2;
                        ylo = f(zp);
                        zp(i) = z(i);
                        hess(i,i) = (yhi + ylo - 2*y)/(eps.^2 / 4);
                        grad(i) = (yhi - ylo) / eps;
                    end
                else
                    grad = this.numerical_jac(f,z);
                    hess = this.numerical_jac(@(z) this.numerical_jac(f,z), z);
                    hess = (hess + hess')/2;
                end
                mineigv = min(eig(hess));
                hess = hess+(mineigv+0.1)*eye(length(z));
            end
        end
        
        % jacobian
        function grad = numerical_jac(this,f,z)
            % modified from Pieter Abbeel's CS287            
            y = f(z);
            
            grad = zeros(length(y), length(z));
            
            eps = 1e-5;
            zp = z;
            
            for i=1:length(z)
                zp(i) = z(i) + eps/2;
                yhi = f(zp);
                zp(i) = z(i) - eps/2;
                ylo = f(zp);
                zp(i) = z(i);
                grad(:,i) = (yhi - ylo) / eps;
            end            
        end
        
         %% sensor model approximation
        % exact value of gamma
        function gam_exact = gam(this,z,theta,x0,alp1,alp2,alp3) 
            gam_exact = 1/this.gam_den(z,theta,x0,alp1,alp2,alp3);
        end
                
        % denominator of gamma
        function gam_d = gam_den(this,z,theta,x0,alp1,alp2,alp3)
            gam_d = this.gam_den1(z,x0,alp1)*this.gam_den2(z,x0,theta,alp2)*this.gam_den3(z,x0,theta,alp3);
        end
        
        function gam_den = gam_den1(this,z,x0,alp)
            tmp = alp*(sum((x0-z).^2)-this.r^2);
            if tmp > this.thr
                gam_den = 1+exp(this.thr);
            else
                gam_den = 1+exp(tmp);
            end
        end
        
        function gam_den = gam_den2(this,z,x0,theta,alp)
            tmp = alp*[sin(theta-this.theta0),-cos(theta-this.theta0)]*(x0-z);
            if tmp > this.thr
                gam_den = 1+exp(this.thr);
            else
                gam_den = 1+exp(tmp);
            end
        end
        
        function gam_den = gam_den3(this,z,x0,theta,alp)
            tmp = alp*[-sin(theta+this.theta0),cos(theta+this.theta0)]*(x0-z);
            if tmp > this.thr
                gam_den = 1+exp(this.thr);
            else
                gam_den = 1+exp(tmp);
            end
        end
        
        % approximate gamma
        function gam_ap = gam_aprx (this,z,theta,x0,z_ref,theta_ref,alp1,alp2,alp3) 
            gam_ap = this.gam(z_ref,theta_ref,x0,alp1,alp2,alp3)...
                +this.gam_grad(z_ref,theta_ref,x0,alp1,alp2,alp3)'*[z-z_ref;theta-theta_ref];
        end
        
        % gradient of gamma
        function gam_g = gam_grad(this,z,theta,x0,alp1,alp2,alp3)
            gam_g = -(this.gam_den1_grad(z,x0,alp1)*this.gam_den2(z,x0,theta,alp2)*this.gam_den3(z,x0,theta,alp3)+...
                this.gam_den1(z,x0,alp1)*this.gam_den2_grad(z,x0,theta,alp2)*this.gam_den3(z,x0,theta,alp3)+...
                this.gam_den1(z,x0,alp1)*this.gam_den2(z,x0,theta,alp2)*this.gam_den3_grad(z,x0,theta,alp3))/...
                (this.gam_den1(z,x0,alp1)*this.gam_den2(z,x0,theta,alp2)*this.gam_den3(z,x0,theta,alp3))^2;
        end
        
        function gam_grad = gam_den1_grad(this,z,x0,alp)
            gam_grad =  [2*(this.gam_den1(z,x0,alp)-1)*alp*(z-x0);0];
        end
        
        function gam_grad = gam_den2_grad(this,z,x0,theta,alp)
            gam_grad =  (this.gam_den2(z,x0,theta,alp)-1)*alp*[-sin(theta-this.theta0);...
                cos(theta-this.theta0);...
                [cos(theta-this.theta0),sin(theta-this.theta0)]*(x0-z)];
        end
        
        function gam_grad = gam_den3_grad(this,z,x0,theta,alp)
            gam_grad =  (this.gam_den3(z,x0,theta,alp)-1)*alp*[sin(theta+this.theta0);...
                -cos(theta+this.theta0);...
                [cos(theta+this.theta0),sin(theta+this.theta0)]*(z-x0)];
        end
        
        % linearized update rule for covariance
        function p = p_aprx(this,z,theta,p1,p2,x0,t1,t2,z_ref,theta_ref,p1_ref,p2_ref,alp1,alp2,alp3) 
            gam_ref = this.gam(z_ref,theta_ref,x0,alp1,alp2,alp3);            
            gam_grad_ref = this.gam_grad(z_ref,theta_ref,x0,alp1,alp2,alp3);
%             display('Robot.m line 1252')
            
%             display(gam_ref)
            if abs(gam_ref) < 10^-4
                gam_ref = 0;
            end
            
%             display(gam_grad_ref)
            tmp_idx = abs(gam_grad_ref) < 10^-4;
            gam_grad_ref(tmp_idx) = 0;
           
            p = p1_ref-gam_ref*(t1*p1_ref+t2*p2_ref)+...
                (1-gam_ref*t1)*(p1-p1_ref)-...
                gam_ref*t2*(p2-p2_ref)-(t1*p1_ref+t2*p2_ref)*...
                gam_grad_ref'*([z-z_ref;theta-theta_ref]);
            
%             p = p1_ref-this.gam(z_ref,theta_ref,x0,alp1,alp2,alp3)*(t1*p1_ref+t2*p2_ref)+...
%                 (1-this.gam(z_ref,theta_ref,x0,alp1,alp2,alp3)*t1)*(p1-p1_ref)-...
%                 this.gam(z_ref,theta_ref,x0,alp1,alp2,alp3)*t2*(p2-p2_ref)-(t1*p1_ref+t2*p2_ref)*...
%                 this.gam_grad(z_ref,theta_ref,x0,alp1,alp2,alp3)'*([z-z_ref;theta-theta_ref]);
        end
        %}
        
        %% %%%%% visualization for debug purpose
        function plotSeedTraj(this,z,x,fld)
            % Plotting for the seed trajectory, i.e. the one
            % fromgenInitState_kf
            
            hdl1 = plot(z(1,:),z(2,:),'k','markers',1);
            set(hdl1,'MarkerFaceColor','k');
            set(hdl1,'MarkerEdgeColor','k');
            set(hdl1,'Color','k');
            set(hdl1,'LineStyle','-.');
%             set(hdl1,'Marker','o');
            
%             hdl2 = plot(x(1,:),x(2,:),'b','markers',1);
%             set(hdl2,'MarkerFaceColor','b');
%             set(hdl2,'MarkerEdgeColor','b');%
%             set(hdl2,'Color','b');
%             %     set(hdl2,'LineStyle','-');
% %             set(hdl2,'Marker','*');
            
            
            % draw FOV
            this.drawFOV(z,fld,'plan')
        end
        
        function plotPlannedTraj(this,z,x,fld)
            % Plotting for the planned trajectory, i.e. the one
            % cvxPlanner_kf
            
            hdl1 = plot(z(1,:),z(2,:),'r','markers',1);
            set(hdl1,'MarkerFaceColor','r');
            set(hdl1,'MarkerEdgeColor','r');
            set(hdl1,'Color','r');
            set(hdl1,'LineStyle','-.');
%             set(hdl1,'Marker','o');
            
%             hdl2 = plot(x(1,:),x(2,:),'b','markers',1);
%             set(hdl2,'MarkerFaceColor','b');
%             set(hdl2,'MarkerEdgeColor','b');%
%             set(hdl2,'Color','b');
%             %     set(hdl2,'LineStyle','-');
% %             set(hdl2,'Marker','*');
            
            
            % draw FOV
            this.drawFOV(z,fld,'plan')
%             for ii = 1:N+1
%                 a1 = z(3,ii)-this.theta0;  % A random direction
%                 a2 = z(3,ii)+this.theta0;
%                 t = linspace(a1,a2,50);
%                 x0 = z(1,ii);
%                 y0 = z(2,ii);
%                 x1 = z(1,ii) + this.r*cos(t);
%                 y1 = z(2,ii) + this.r*sin(t);
%                 plot([x0,x1,x0],[y0,y1,y0],'r--','LineWidth',0.5)
%             end
%             
%             xlim([fld.fld_cor(1),fld.fld_cor(2)]);
%             ylim([fld.fld_cor(3),fld.fld_cor(4)]);
%             box on
%             axis equal
%             drawnow
        end
        
        function drawFOV(this,z,fld,fov_mode)
            if strcmp(fov_mode,'plan')
                % draw the FOV for planned robot state
                N = this.mpc_hor;
                for ii = 1:N+1
                    a1 = z(3,ii)-this.theta0;  % A random direction
                    a2 = z(3,ii)+this.theta0;
                    t = linspace(a1,a2,50);
                    x0 = z(1,ii);
                    y0 = z(2,ii);
                    x1 = z(1,ii) + this.r*cos(t);
                    y1 = z(2,ii) + this.r*sin(t);
                    plot([x0,x1,x0],[y0,y1,y0],'r--','LineWidth',0.5)
                end
            elseif strcmp(fov_mode,'cur')
                % draw the FOV for current robot state
                a1 = z(3)-this.theta0;  % A random direction
                a2 = z(3)+this.theta0;
                t = linspace(a1,a2,50);
                x0 = z(1);
                y0 = z(2);
                x1 = z(1) + this.r*cos(t);
                y1 = z(2) + this.r*sin(t);
                plot([x0,x1,x0],[y0,y1,y0],'r--','LineWidth',0.5)
            end
            xlim([fld.fld_cor(1),fld.fld_cor(2)]);
            ylim([fld.fld_cor(3),fld.fld_cor(4)]);
            box on
            axis equal
            drawnow
        end
    end
end