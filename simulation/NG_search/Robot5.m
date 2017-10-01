 classdef Robot5
     % this class uses ipopt for solution
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
        snum;
        
        
        % performance metrics
        ml_pos;
        ent_pos;
    end
    
    methods
        function this = Robot5(inPara)
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
            % opt 1: plain random resampling
            %{
            idx = randsample(1:np, np, true, w);
            new_particles = pred_par(:,idx);
            this.particles = new_particles;
            %}
            
            % opt 2: low variance sampling
            %
            M = 1/np;
            U = rand(1)*M;
            new_particles = zeros(size(pred_par));
            tmp_w = w(1);
            ii = 1;
            jj = 1;
            while (jj <= np)
                while (tmp_w < U+(jj-1)*M)
                    ii = ii+1;
                    tmp_w = tmp_w+w(ii);                    
                end
                new_particles(:,jj) = pred_par(:,ii);
                jj = jj + 1;
            end
            this.particles = new_particles;
            %}
            
            %% gmm fitting
            max_gmm_num = this.max_gmm_num; % maximum gmm component number
            
            gmm_model = cell(max_gmm_num,1);
            opt = statset('MaxIter',1000);
            AIC = zeros(max_gmm_num,1);
%             tic;
            for kk = 1:max_gmm_num
                gmm_model{kk} = fitgmdist(new_particles',kk,'Options',opt,...
                    'Regularize',0.001,'CovType','full');
                AIC(kk)= gmm_model{kk}.AIC;
            end
%             fprintf('Robot.m, Line %d. gmm fitting takes time as:',MFileLineNr());
%             toc;
            
            [minAIC,numComponents] = min(AIC);
            
            best_model = gmm_model{numComponents};
            this.gmm_num = numComponents;
            this.gmm_mu = best_model.mu';
            this.gmm_sigma = best_model.Sigma;
            % optimization will generate bad behavior if sigma is too
            % small. so constrain the min sig
            for ii = 1:this.gmm_num
                mineigv = min(eig(this.gmm_sigma(:,:,ii)));
                if mineigv < 1
                    this.gmm_sigma(:,:,ii) = 1*eye(2);
                end
            end
            
            this.wt = best_model.PComponents';
           
            
            % ad-hoc way of merging close gmm componenets
            % needs a more systematic way of doing this later
            %{
            for kk = 1:this.gmm_num
                if kk == this.gmm_num
                    break
                end
                [idx,cetr,sumd] = kmeans(this.gmm_mu',kk);
                if all(sumd < 2)
                    break
                end
            end
            
            if kk < this.gmm_num
               new_mu = zeros(2,kk);
               new_sigma = zeros(2,2,kk);
               new_wt = zeros(kk,1);
               tmp_P = reshape(this.gmm_sigma,4,this.gmm_num);
               for ii = 1:kk
                   new_mu(:,ii) = sum(this.wt(idx==ii)'.*this.gmm_mu(:,idx==ii),2)/sum(this.wt(idx==ii));                   
                   new_sigma(:,:,ii) = reshape(sum(this.wt(idx==ii)'.*tmp_P(:,idx==ii),2)/sum(this.wt(idx==ii)),2,2);
                   new_wt(ii) = sum(this.wt(idx==ii));
               end
               this.gmm_mu = new_mu;
               this.gmm_sigma = new_sigma;
               this.wt = new_wt;
               this.gmm_num = kk;
            end
            %}
            
            % convert the data form to be compatible with the main code
            this.est_pos = this.gmm_mu(:);
            for ii = 1:this.gmm_num %numComponents
                this.P{ii} = this.gmm_sigma(:,:,ii);
            end
        end
            
       
        %% planning
        
        % try re-writing the problem using cvx solver. Different from
        % cvxPlanner below, which formulates the problem as a convex P (turns
        % out not!), this one formulates the problem as QP each iteration
        function [optz,optu, s, snum, merit, model_merit, new_merit] = cvxPlanner_ipopt(this,fld,optz,optu,plan_mode) % cvxPlanner(this,fld,init_sol)            
            % use ipopt to solve inner layer optimization. The outter layer
            % changes alpha for approximating gamma
            
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
            Pinverse_ref = init_sol.Pinverse;       
            auxt_ref = init_sol.auxt;
            auxm_ref = init_sol.auxm;
            gam_var = init_sol.gam_var;
            slkm = init_sol.slkm;
            slkt = init_sol.slkt;
            auxgau_ref = init_sol.auxgau;
            
            % construct the state
            % s is the collection of all decision variables
            % s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:)]   
            [s,snum] = this.setState(zref,uref,xref,x_pred_ref,Pref,P_pred_ref,Kref,Pinverse_ref,auxt_ref,auxm_ref,gam_var,slkm,slkt,auxgau_ref);
            
%             cfg.snum = snum;
%             this.cfg = cfg;
%             this.snum = snum;
            
            %%% functions and parameters for scp                     
            %%% objective
%             obj = @(s) this.getObj(s,snum);
            
            % Q and q in quad linear obj term: x'*Q*x/2+q'*x
            fQ = zeros(length(s)); 
            q = zeros(1,length(s));
            
            %%% kinematics constraints
            objKin = @(s) this.getKinConstr(s,snum);
            
            %%% constraint for P and Pinverse
            objPinverse = @(s) this.getPinverseConstr(s,snum);
            
            %%% constraint for auxiliary variables
            objAux = @(s) this.getAuxConstr(s,snum);
            
            %%% linear equality constraints          
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
            constrLinIneq = @(s) this.getLinIneqConstr(s,snum,fld);
            
%             penalty_coeff = cfg.initial_penalty_coeff; % Coefficient of l1 penalties 
%             trust_box_size = cfg.initial_trust_box_size; % The trust region will be a box around the current iterate x.            
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
                        objBel = @(s) this.getBelConstr(s,snum,alp1,alp2,alp3); % belief update
                end
                
                %%% objective
                obj = @(s) this.getObj(s,snum,alp1,alp2,alp3); 
                
                %% layer 2: rely on ipopt itself
                while (1)
                    objNLineq = @(s) 0; % objAux(s); % nonlinear inequality constraint here
                    objNLeq = @(s) [objKin(s);objBel(s);objPinverse(s);objAux(s)]; %% nonlinear equality constraint here
                    % 
                    % use yalmip to define
                    sp = sdpvar(length(s),1);
                    % define auxiliary variables to simplify the obj for yalmip
                    
                    constr = [constrLinIneq(sp) <= 0; constrLinEq(sp) == 0;...
                        objNLineq(sp) <= 0; objNLeq(sp) == 0];
                    assign(sp,s);
                    opt = sdpsettings('verbose',5,'solver','ipopt','usex0',1,'debug',1);
                    opt.ipopt.tol = 10^-3;
%                     opt.ipopt.hessian_approximation = 'exact';
%                     opt.ipopt.mu_strategy = 'monotone';
                    
                    optobj = obj(sp);
                    
                    sol = optimize(constr,optobj,opt);
                    
                    s = value(sp);
                    
                    if sol.problem == 0
                        success = true;
                        %                     s = value(sp);
                        break
                    else
                        success = false;
                        fprintf('  [ipopt loop] Robot.m, line %d.\n',MFileLineNr())
                        fprintf('  solution failed\n')
                        fprintf('  the issue is %d\n: %s\n',sol.problem,sol.info)
%                         break
                    end
                end
                %%% layer 2 ends
                
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
                
                % update initial solution for next iteration
                s = this.updOptVar(s,snum,alp1,alp2,alp3);
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                if max(tmp_dif) <= cfg.gamma_tol %0.05
                    fprintf('  [gamma loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('  Gamma error is within tolerance, break out of gamma loop\n')
                    break
                elseif gam_iter >= cfg.max_gam_iter 
                    cprintf('Red',sprintf('  [gamma loop] Robot.m, line %d.  max gamma iteration reached. break out of gamma loop\n',MFileLineNr()))
                    break
                else
                    cprintf('Magenta',sprintf('  [gamma loop] Robot.m, line %d.  gamma_exact is not close, go to another iteration.\n',MFileLineNr()))
                    alp1 = alp1*alp_inc;
                    alp2 = alp2*alp_inc;
                    alp3 = alp3*alp_inc;
                    % (added on 9/26) I think s should be updated if alp's
                    % are updated
                    s = this.updOptVar(s,snum,alp1,alp2,alp3);
                end
                gam_iter = gam_iter+1;
            end % loop 1 ends
            
            uref = this.convState(s,snum,'u');
            
            % use actual dynamics to simulate
            zref = this.simState(uref,zref(:,1));
            
            %%% also revise this part in sqp code
%             if ~success %if not solved successfully, reuse previous input
%                 uref = [uref(:,2:end),uref(:,end)];
%             end
            
            optz = zref;
            optu = uref;
            
            % visualize the FOV along the planned path
            %%% xref in this part needs change when infeasibility happens
            this.plotPlannedTraj(optz,xref,fld)
%             }

            merit = 0;
            model_merit = 0; 
            new_merit = 0;
        end
        
        function [optz,optu, s, snum, merit, model_merit, new_merit] = cvxPlanner_scp(this,fld,optz,optu,plan_mode) % cvxPlanner(this,fld,init_sol)
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
            Pinverse_ref = init_sol.Pinverse;       
            auxt_ref = init_sol.auxt;
            auxm_ref = init_sol.auxm;
            gam_var = init_sol.gam_var;
            slkm = init_sol.slkm;
            slkt = init_sol.slkt;
            auxgau_ref = init_sol.auxgau;
            
            % construct the state
            % s is the collection of all decision variables
            % s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:)]   
            [s,snum] = this.setState(zref,uref,xref,x_pred_ref,Pref,P_pred_ref,Kref,Pinverse_ref,auxt_ref,auxm_ref,gam_var,slkm,slkt,auxgau_ref);
            
%             cfg.snum = snum;
%             this.cfg = cfg;
%             this.snum = snum;
            
            %%% functions and parameters for scp                     
            %%% objective
%             obj = @(s) this.getObj(s,snum);
            
            % Q and q in quad linear obj term: x'*Q*x/2+q'*x
            fQ = zeros(length(s)); 
            q = zeros(1,length(s));
            
            %%% kinematics constraints
            objKin = @(s) this.getKinConstr(s,snum);
            
            %%% constraint for P and Pinverse
            objPinverse = @(s) this.getPinverseConstr(s,snum);
            
            %%% constraint for auxiliary variables
            objAux = @(s) this.getAuxConstr(s,snum);
            
            %%% linear equality constraints          
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
            constrLinIneq = @(s) this.getLinIneqConstr(s,snum,fld);
            
%             penalty_coeff = cfg.initial_penalty_coeff; % Coefficient of l1 penalties 
%             trust_box_size = cfg.initial_trust_box_size; % The trust region will be a box around the current iterate x.            
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
                        objBel = @(s) this.getBelConstr(s,snum,alp1,alp2,alp3); % belief update
                end
                
                %%% objective
                obj = @(s) this.getObj(s,snum,alp1,alp2,alp3); 
                
                %% layer 2: penalty iteration
                penalty_coeff = cfg.initial_penalty_coeff; % Coefficient of l1 penalties  
                
                %                 ctol = 1e-2;
                %                 mu = 0.1;
                pen_iter = 0;
                while(1) % loop 2 starts
                    pen_iter = pen_iter+1;
                    fprintf('    [Penalty loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('    pentaly loop #%d\n',pen_iter);
                    
                    trust_box_size = cfg.initial_trust_box_size; % The trust region will be a box around the current iterate x.
                    
                    objNLineq = @(s) 0; % nonlinear inequality constraint here
                    objNLeq = @(s) [objKin(s);objBel(s);objPinverse(s);objAux(s)]; % nonlinear equality constraint here
                    %                     objNLeq = @(s) 0;
                    
                    %% loop 3: trust region SQP
                    A_ineq = [];
                    b_ineq = [];
                    A_eq = [];
                    b_eq = [];
                    
                    % make initial solution feasible
                    [s,success] = this.find_closest_feasible_point(s,constrLinIneq,constrLinEq);
                    if (~success)
                        return;
                    end
                    
                    % loop 3 starts
                    [s, trust_box_size, success, merit, model_merit, new_merit] = this.minimize_merit_function(s, fQ, q, obj, A_ineq, b_ineq, A_eq, b_eq, constrLinIneq, constrLinEq, objNLineq, objNLeq, hinge, abssum, cfg, penalty_coeff, trust_box_size, snum, fld);
                    % loop 3 ends
                    
                    if(hinge(objNLineq(s)) + abssum(objNLeq(s)) < cfg.cnt_tolerance || pen_iter >= cfg.max_penalty_iter ) %cfg.max_iter
                        break;
                    end
                    %                     trust_box_size = cfg.initial_trust_box_size;
                    penalty_coeff = cfg.merit_coeff_increase_ratio*penalty_coeff;
                end % loop 2 ends
                
                if ~success %infea_flag
                    % if CVS infeasible, directly reuse solution from
                    % previous step
                    break
                end
                %%% layer2: sqp solver ends
                
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
                
                % update initial solution for next iteration
                s = this.updOptVar(s,snum,alp1,alp2,alp3);
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
                if max(tmp_dif) <= cfg.gamma_tol %0.05
                    fprintf('  [gamma loop] Robot.m, line %d\n', MFileLineNr())
                    fprintf('  Gamma error is within tolerance, break out of gamma loop\n')
                    break
                elseif gam_iter >= cfg.max_gam_iter 
                    cprintf('Red',sprintf('  [gamma loop] Robot.m, line %d.  max gamma iteration reached. break out of gamma loop\n',MFileLineNr()))
                    break
                else
                    cprintf('Magenta',sprintf('  [gamma loop] Robot.m, line %d.  gamma_exact is not close, go to another iteration.\n',MFileLineNr()))
                    alp1 = alp1*alp_inc;
                    alp2 = alp2*alp_inc;
                    alp3 = alp3*alp_inc;
                    % (added on 9/26) I think s should be updated if alp's
                    % are updated
                    s = this.updOptVar(s,snum,alp1,alp2,alp3);
                end
                gam_iter = gam_iter+1;
            end % loop 1 ends
            
            uref = this.convState(s,snum,'u');
            
            % use actual dynamics to simulate
            zref = this.simState(uref,zref(:,1));
            
            %%% also revise this part in sqp code
%             if ~success %if not solved successfully, reuse previous input
%                 uref = [uref(:,2:end),uref(:,end)];
%             end
            
            optz = zref;
            optu = uref;
            
            % visualize the FOV along the planned path
            %%% xref in this part needs change when infeasibility happens
            this.plotPlannedTraj(optz,xref,fld)
%             }

            merit = 0;
            model_merit = 0; 
            new_merit = 0;
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
            
            alp1 = this.alp1;
            alp2 = this.alp2;
            alp3 = this.alp3;
            
            % open-loop prediction of target position. no measurement info 
            % is used
            % get x, x_pred
            x_olp = zeros(2*this.gmm_num,N+1);
            x_olp(:,1) = this.est_pos(:);
            for ii = 1:N
                for jj = 1:this.gmm_num
                    x_olp(2*jj-1:2*jj,ii+1) = f(x_olp(2*jj-1:2*jj,ii));
                end
            end            
            init_state.x_olp = x_olp;
            init_state.x = x_olp;
            init_state.x_pred = x_olp(:,2:end);                        
            
            % get z,u
            if isempty(prev_state)
                % if this is the 1st time of calling ngPlanner, use motion
                % primitives
                
                u = [0;0.1]*ones(1,N);
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
            
%             % get theta_bar
%             theta_bar = zeros(this.gmm_num,N+1);
%             for ii = 1:N+1
%                 for jj = 1:this.gmm_num
%                     tmp_vec = x_olp(2*jj-1:2*jj,ii)-z(1:2,ii);
%                     theta_bar(jj,ii) = atan2(tmp_vec(2),tmp_vec(1));
%                 end
%             end
%             init_state.theta_bar = theta_bar;
            
            % get P, P_pred, K
            P = zeros(2,2,this.gmm_num,N+1);
            P_pred = zeros(2,2,this.gmm_num,N);
            K = zeros(2,2,this.gmm_num,N);
            Pinverse = zeros(size(P)); % used for ipopt since yalmip does not allow inverse of matrix variables
            gam_var = zeros(this.gmm_num,N);
            
            for jj = 1:this.gmm_num
                P(:,:,jj,1) = this.P{jj};
            end
            
            for ii = 1:N                
                for jj = 1:this.gmm_num
                    A = del_f(x_olp(2*jj-1:2*jj,ii+1));
                    C = del_h(x_olp(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    P_pred(:,:,jj,ii) = A*P(:,:,jj,ii)*A'+Q;
                    T = C*P_pred(:,:,jj,ii)*C'+R;
                    K(:,:,jj,ii)= P_pred(:,:,jj,ii)*C'/T;
%                     gam = this.inFOV(x_olp(2*jj-1:2*jj,ii));
                    tmp_sum =zeros(2,2);
                    for ll = jj %1:this.gmm_num
                        gam = this.gam_eps(z(1:2,ii+1),z(3,ii+1),x_olp(2*ll-1:2*ll,ii+1),alp1,alp2,alp3);
                        gam_var(ll,ii) = gam;
                        tmp_sum = tmp_sum+this.wt(ll)*gam*K(:,:,jj,ii)*C*P_pred(:,:,jj,ii);
                    end
                    P(:,:,jj,ii+1) = P_pred(:,:,jj,ii)-tmp_sum;
                end
            end
            
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    Pinverse(:,:,jj,ii) = eye(2)/P(:,:,jj,ii);
                end
            end
            
            init_state.P = P;
            init_state.P_pred = P_pred;
            init_state.K = K;    
            init_state.Pinverse = Pinverse; 
            init_state.gam_var = gam_var;  
            
            % auxiliary variable
            t = zeros(this.gmm_num,N);
            m = zeros(this.gmm_num,this.gmm_num,N);
            gau = zeros(this.gmm_num,this.gmm_num,N);
            
            % m
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
                        m(jj,ll,ii-1) = this.wt(ll)*sqrt(det(Pinverse(:,:,ll,ii)))/(2*pi)*exp(-(x_olp(2*jj-1:2*jj,ii)-x_olp(2*ll-1:2*ll,ii))'*Pinverse(:,:,ll,ii)*(x_olp(2*jj-1:2*jj,ii)-x_olp(2*ll-1:2*ll,ii))/2);
                        gau(jj,ll,ii-1) = (x_olp(2*jj-1:2*jj,ii)-x_olp(2*ll-1:2*ll,ii))'*Pinverse(:,:,ll,ii)*(x_olp(2*jj-1:2*jj,ii)-x_olp(2*ll-1:2*ll,ii))/2;
                    end
                end
            end
            
            % t
            for ii = 1:N
                for jj = 1:this.gmm_num
                    t(jj,ii) = -this.wt(jj)*log(sum(m(jj,:,ii)));
                end
            end
            
            init_state.auxm = m;
            init_state.auxt = t;
            init_state.auxgau = gau;
            
            % slack variables
            slkm = zeros(this.gmm_num,this.gmm_num,N);
            slkt = zeros(this.gmm_num,N);
            
            init_state.slkm = slkm;
            init_state.slkt = slkt;
            
            % visualize the seed path
            hold on
%             this.plotSeedTraj(z,x_olp,fld)
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
            % use u to compute state z, which replaces current z to be the
            % new state
            st = this.state;
            this.optu = [this.optu,u(:,1)];
            dt = this.dt;
            this.state = st+[st(4)*cos(st(3));st(4)*sin(st(3));u(:,1)]*dt;
            this.traj = [this.traj,this.state];
        end               
        
        function z = simState(this,u,z0)
            % use u to compute state z. 1st element in z is the current
            % state
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
                    [fgrad, ~] = this.getNumGradHess(f,x,cfg.full_hessian);
%                     [fgrad, fhess] = this.getNumGradHess(f,x,cfg.full_hessian);
                    % diagonal adjustment
%                     mineig = min(eigs(fhess));
%                     if mineig < 0
%                         fprintf('[sqp loop] Robot.m, line %d', MFileLineNr())
%                         fprintf('negative hessian detected. adjusting by %.3g\n',-mineig);
%                         fhess = fhess + eye(dim_x) * ( - mineig);
%                     end
                else %if ~isempty(f)
%                     [fval, fgrad, fhess] = f(x);
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
                    
                    
                    cvx_begin quiet
                        variables xp(dim_x, 1);
                        %                 minimize fquadlin(x) + fval + (x'*Q'+q+fgrad)*(xp-x) + penalty_coeff*( hinge(gval+gjac*(xp-x)) + abssum(hval+hjac*(xp-x)) )
                        minimize fquadlin(xp) + fval + (fgrad)*(xp-x) + penalty_coeff*( hinge(gval+gjac*(xp-x)) + abssum(hval+hjac*(xp-x)) )
                        subject to
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
                fprintf('    Robot.m, line %d\n', MFileLineNr())
                fprintf('    initialization doesn''t satisfy linear constraints. finding the closest feasible point\n');
                
                cvx_begin quiet
                variables('x(length(x0))')
                minimize('sum((x-x0).^2)')
                subject to
                    hlin(x) == 0;
                    glin(x) <= 0;
                cvx_end
                
                if strcmp(cvx_status,'Failed') || strcmp(cvx_status,'Infeasible')
                    success = false;
                    fprintf('    Robot.m, line %d\n', MFileLineNr())
                    fprintf('    Couldn''t find a point satisfying linear constraints\n')
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
            x_orig = this.convState(s,snum,'x'); 
            xpred_orig = this.convState(s,snum,'x_pred');
            P_orig = this.convState(s,snum,'P');
            P_pred_orig = this.convState(s,snum,'P_pred');
            K_orig = this.convState(s,snum,'K');
            gam_var_orig = this.convState(s,snum,'gam_var');
            slkm = this.convState(s,snum,'slkm');
            slkt = this.convState(s,snum,'slkt');
            
            N = this.mpc_hor;
            
%             alp1 = this.alp1;
%             alp2 = this.alp2;
%             alp3 = this.alp3;
            
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            f = this.target.f;
            del_f = this.target.del_f;
            Q = this.target.Q;

            % z
            z = this.simState(u_orig,z_orig(:,1));
            
            % 9/6: this commented part is used for the KF case. A new code 
            % that is intended to be general is written. But I need to test
            % if the new part can work for KF or not. If it can, then
            % remove this commented part.
            %{
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
            %}
            
            % new, generat code
            % belief state
            x = zeros(size(x_orig));
            P = zeros(size(P_orig));
            K = zeros(size(K_orig));
            xpred = zeros(size(xpred_orig));
            P_pred = zeros(size(P_pred_orig));
            Pinverse = zeros(size(P_orig));
            gam_var = zeros(size(gam_var_orig));
            
            % initialization
            x(:,1) = x_orig(:,1);
            %%%%%% fix this in sqp version
            for jj = 1:this.gmm_num
                P(:,:,jj,1) = P_orig(:,:,jj,1);
            end
            
%             for ii = 1:N
%                 for jj = 1:this.gmm_num
%                     x(2*jj-1:2*jj,ii+1) = f(x(2*jj-1:2*jj,ii));
%                 end
%             end
            
            for ii = 1:N
                % predicte all variables
                for jj = 1:this.gmm_num
                    A = del_f(x(2*jj-1:2*jj,ii+1));                    
                    xpred(2*jj-1:2*jj,ii) = f(x(2*jj-1:2*jj,ii));
                    P_pred(:,:,jj,ii) = A*P(:,:,jj,ii)*A'+Q;                    
                end
                
                % compute variables using measurement
                for jj = 1:this.gmm_num
                    C = del_h(x(2*jj-1:2*jj,ii+1),z(1:2,ii+1));
                    T = C*P_pred(:,:,jj,ii)*C'+R;
                    K(:,:,jj,ii)= P_pred(:,:,jj,ii)*C'/T;
                    %                 gam = this.inFOV(x_olp(:,ii));
                                        
                    gm = zeros(this.gmm_num,1);
                    tmp = 0;
                    for ll = 1:this.gmm_num
                        gm(ll) = this.gam(z(1:2,ii+1),z(3,ii+1),xpred(2*ll-1:2*ll,ii),alp1,alp2,alp3);
                        tmp = tmp+this.wt(ll)*gm(ll)*K(:,:,jj,ii)*...
                            (h(xpred(2*ll-1:2*ll,ii),z(1:2,ii+1))-h(xpred(2*jj-1:2*jj,ii),z(1:2,ii+1)));
                        gam_var(ll,ii) = gm(ll);
                        
                    end
                    x(2*jj-1:2*jj,ii+1) = xpred(2*jj-1:2*jj,ii)+tmp;
                    P(:,:,jj,ii+1) = P_pred(:,:,jj,ii)-this.wt'*gm*K(:,:,jj,ii)*C*P_pred(:,:,jj,ii);                   
                end
            end
            
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    Pinverse(:,:,jj,ii) = eye(2)/P(:,:,jj,ii);
                end
            end
            
            % auxiliary variables
            auxt = zeros(this.gmm_num,N);
            auxm = zeros(this.gmm_num,this.gmm_num,N);
            auxgau = zeros(this.gmm_num,this.gmm_num,N);
            
            % auxm
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
                        auxm(jj,ll,ii-1) = this.wt(ll)*sqrt(det(Pinverse(:,:,ll,ii)))/(2*pi)*exp(-(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))/2);
                        auxgau(jj,ll,ii-1) = (x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))/2;
                    end
                end
            end
            
            % auxt
            for ii = 1:N
                for jj = 1:this.gmm_num
                    auxt(jj,ii) = -this.wt(jj)*log(sum(auxm(jj,:,ii)));
                end
            end
            
            snew = setState(this,z,u_orig,x,xpred,P,P_pred,K,Pinverse,auxt,auxm,gam_var,slkm,slkt,auxgau);
        end
        
        %% objective function
        % for general NGP. 
        function val = getObj(this,s,snum,alp1,alp2,alp3)
            x = this.convState(s,snum,'x');
            u = this.convState(s,snum,'u');
            P = this.convState(s,snum,'P');
            z = this.convState(s,snum,'z');
            Pinverse = this.convState(s,snum,'Pinverse');
            auxt = this.convState(s,snum,'auxt');
            slkm = this.convState(s,snum,'slkm'); 
            slkt = this.convState(s,snum,'slkt'); 
            
            N = this.mpc_hor;
            val = 0;
            [~,max_idx] = max(this.wt);
            for ii = 2:N+1
               %{               
               for jj = 1:this.gmm_num
                   
                   tmp = 0;
                   for ll = 1:this.gmm_num                       
%                        P(:,:,ll,ii) = (P(:,:,ll,ii)+P(:,:,ll,ii)')/2; % numerical issues happen that makes P non symmetric
                      
                       %%% (9/20) note: this commented out snippet coudl be useful
                       %%% to ensure the robustness of the solver. but
                       %%% could give bad result if P has negative evalue
                       %%% and then compensated to be positive. Leave this
                       %%% snippet for record purpose only
%                        mineigv = min(eig(P(:,:,ll,ii)));                       
%                        if mineigv <= 0
%                            P(:,:,ll,ii) = P(:,:,ll,ii)+(-mineigv+0.01)*eye(size(P(:,:,ll,ii),1));
%                        end                       
%                        tmp = tmp+this.wt(ll)*mvnpdf(x(2*jj-1:2*jj,ii),x(2*ll-1:2*ll,ii),P(:,:,ll,ii));
%                         tmp = tmp+this.wt(ll)*exp(-(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))/2)...
%                            *sqrt(Pinverse(1,1,ll,ii)*Pinverse(2,2,ll,ii)-Pinverse(1,2,ll,ii)*Pinverse(2,1,ll,ii))/(2*pi);
                        
                   end
                   tmp_dis = sum((x(2*jj-1:2*jj,ii)-z(1:2,ii)).^2);
                   tmp_dis = abs(sum((x(2*jj-1:2*jj,ii)-z(1:2,ii)).^2)-1); % distance between sensor and MLE target position
                   
                   % added on 8/20
                   % temporarily add the diff between gamma and 1 as
                   % penalty
%                    tar_pos = x(2*jj-1:2*jj,ii);                  

                   val = val-10*this.wt(jj)*log(tmp);%+this.wt(jj)*tmp_dis; %+(1-this.gam(z(1:2,ii),z(3,ii),tar_pos,alp1,alp2,alp3));                   
               end 
               %}
               val = val+sum(auxt(:,ii-1));   
               
               val = val+sum(u(:,ii-1).^2); % penalize on control input
               val = val+sum((x(2*max_idx-1:2*max_idx,ii)-z(1:2,ii)).^2);
%                val = val+abs(sum((x(2*max_idx-1:2*max_idx,ii)-z(1:2,ii)).^2)-2); % penalize the distance between sensor and MLE target postion with maximal weight
            end
            val = val+sum(slkm(:).^2)+sum(slkt(:).^2);
            
%             val = 0.1*val;
%             val = 0;
        end    

        function val = getQuadObj(this,x,snum)
            % linear quadratic term in the objective function
            val = 0;
%             u = this.convState(x,snum,'u');
%             val = sum(sum(u.^2));
        end
    
        %% constraints
        % for general NGP
        % NL equality constraints
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
        
        function h = getBelConstr(this,s,snum,alp1,alp2,alp3)
            % belief update constraint (nonlinear equality)
            z = this.convState(s,snum,'z'); 
            u = this.convState(s,snum,'u'); 
            x = this.convState(s,snum,'x'); 
            xpred = this.convState(s,snum,'x_pred'); 
            P = this.convState(s,snum,'P'); 
            P_pred = this.convState(s,snum,'P_pred');
            K = this.convState(s,snum,'K');
            gam_var = this.convState(s,snum,'gam_var');
%             h = 0;
            
            val1 = [];
            val2 = [];
            val3 = [];
            val4 = [];
            val5 = [];
            
            N = this.mpc_hor;            
            dt = this.dt;
            gam = @this.gam;
            gam_den = @this.gam_den;
            h = this.h;
            del_h = this.del_h;
            R = this.R;
            
%             alp1 = this.alp1;
%             alp2 = this.alp2;
%             alp3 = this.alp3;
            
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
                    tmp = (K(:,:,jj,ii)*(C(:,:,jj,ii)*P_pred(:,:,jj,ii)*C(:,:,jj,ii)'+R)-...
                        P_pred(:,:,jj,ii)*C(:,:,jj,ii)'); % used a scaling factor 
                    val1 = [val1;tmp(:)];
                    
                    % mean
                    %%%%% note: this part is not used currently, since we
                    %%%%% assume MAP. In future version, I will probably
                    %%%%% use this part once MAP assumption is removed
                    %{
                    % x_k+1|k+1=x_k+1|k+\sum
                    % w_ll*gamma_ll*K*(h(x^ll_k+1|k)-h(x_k+1|k))
                    gm = sdpvar(this.gmm_num,1);
                    
                    tmp2 = 0;
                    for ll = 1:this.gmm_num
                        gm(ll) = this.gam(z(1:2,ii+1),z(3,ii+1),xpred(2*ll-1:2*ll,ii),alp1,alp2,alp3);
                        tmp2 = tmp2+this.wt(ll)*gm(ll)*K(:,:,jj,ii)*...
                            (h(xpred(2*ll-1:2*ll,ii),z(1:2,ii+1))-h(xpred(2*jj-1:2*jj,ii),z(1:2,ii+1)));
                    end
                    val3 = [val3;x(2*jj-1:2*jj,ii+1)-xpred(2*jj-1:2*jj,ii)-tmp2];
                    %}
                    
                    % covariance
                    % P_k+1|k+1=P_k+1|k-\sum w_ll*gamma_ll*K*C*P_k+1|k
                    tmp_sum = zeros(2,2);
                    T = K(:,:,jj,ii)*C(:,:,jj,ii);
                    
                    for ll = jj%1:this.gmm_num
                        %%% (9/27) for now I just assume ll==jj.
                        %%% note: gamma depends on ll, C depends
                        %%% only on jj
%                         tar_pos = xpred(2*ll-1:2*ll,ii); % use each gmm component mean as a possible target position
%                         tmp_sum = tmp_sum+this.wt(ll)*gam(z(1:2,ii+1),z(3,ii+1),...
%                         tar_pos,alp1,alp2,alp3)*T*P_pred(:,:,jj,ii);
                        
                        tmp_sum = tmp_sum+this.wt(ll)*gam_var(ll,ii)*T*P_pred(:,:,jj,ii);                                                
                    end                                        
%                     h = h+sum(sum(abs(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum)));
                    tmp2 = triu(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)+tmp_sum);
                    val2 = [val2;tmp2(:)];
                    
                    % constraint for gam_var
                    % gam_var=1/gam_den
                    % the reason that gam_var is defined outside ll for
                    % loop is because when ll=1:this.gmm_num, gam_var will
                    % be repeated constrained, resulting in degenerate
                    % constraints
                    
%                         val4 = [val4;gam_var(ll,ii)*gam_den(z(1:2,ii+1),z(3,ii+1),...
%                             tar_pos,alp1,alp2,alp3)-1];
%                         val4 = [val4;log(gam_var(ll,ii))+log(gam_den(z(1:2,ii+1),z(3,ii+1),...
%                             tar_pos,alp1,alp2,alp3))];
%                         val4 = [val4;log(gam_var(ll,ii))+log(this.gam_den1(z(1:2,ii+1),...
%                             tar_pos,alp1))+log(this.gam_den2(z(1:2,ii+1),...
%                             tar_pos,z(3,ii+1),alp2))+log(this.gam_den3(z(1:2,ii+1),...
%                             tar_pos,z(3,ii+1),alp3))];

                    % gam_var*\prod_i(1+exp^powi)=\prod_i(eps+exp^powi)                    
                    eps1 = 10^-3;
                    tar_pos = xpred(2*jj-1:2*jj,ii); % use each gmm component mean as a possible target position
                    [~,pow1] = this.gam_den1(z(1:2,ii+1),tar_pos,alp1);
                    [~,pow2] = this.gam_den2(z(1:2,ii+1),tar_pos,z(3,ii+1),alp2);
                    [~,pow3] = this.gam_den3(z(1:2,ii+1),tar_pos,z(3,ii+1),alp3);
                    tmp4 = log(gam_var(jj,ii))+log(1+exp(-pow1))...
                        +log(1+exp(-pow2))+log(1+exp(-pow3))-(log(eps1+exp(-pow1))...
                        +log(eps1+exp(-pow2))+log(eps1+exp(-pow3)));
                    val4 = [val4;tmp4];
                    
                end
            end
            
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    val5 = [val5;-P(1,1,jj,ii)*P(2,2,jj,ii)+...
                        P(1,2,jj,ii)*P(2,1,jj,ii)];
                end
            end
            
            h = [val1;val2;val4;val5];
        end
        
        function h = getPinverseConstr(this,s,snum)
            % enforce that Pinverse is the inverse of P
            N = this.mpc_hor;
            
            P = this.convState(s,snum,'P'); 
            Pinverse = this.convState(s,snum,'Pinverse'); 
            
            h = [];
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    tmp = triu(P(:,:,jj,ii)*Pinverse(:,:,jj,ii))-eye(2);
                    h = [h;tmp(:)];
                end
            end
        end
        
        % NL inequality constraints
        function h = getAuxConstr(this,s,snum)
            x = this.convState(s,snum,'x');
            P = this.convState(s,snum,'P'); 
            Pinverse = this.convState(s,snum,'Pinverse'); 
            auxt = this.convState(s,snum,'auxt'); 
            auxm = this.convState(s,snum,'auxm'); 
            slkm = this.convState(s,snum,'slkm'); 
            slkt = this.convState(s,snum,'slkt'); 
            auxgau = this.convState(s,snum,'auxgau'); 
            
            N = this.mpc_hor;
            
            val1 = [];
            val2 = [];
            val3 = [];
            
            %% (x(jj)-x(ll))'*Pinv(ll)*(x(jj)-x(ll))+2log(2pi/wt(ll))+log(|P(ll)|)+2log(auxm(jj,ll)) <= 0
            % log form
            %{
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
                        tmp = (x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*...
                            Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))+...
                            2*log(2*pi/this.wt(ll))+...
                            log(P(1,1,ll,ii)*P(2,2,ll,ii)-P(1,2,ll,ii)*P(2,1,ll,ii))+...
                            2*log(auxm(jj,ll,ii-1))+slkm(jj,ll,ii-1);
                        val1 = [val1;tmp];
                    end                    
                end
            end
            %}
            
            % exp form
            %{
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
                        tmp = -this.wt(ll)/2*pi*sqrt(Pinverse(1,1,ll,ii)*Pinverse(2,2,ll,ii)...
                            -Pinverse(1,2,ll,ii)*Pinverse(2,1,ll,ii))*...
                            exp(-(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*...
                            Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))/2)+...
                            auxm(jj,ll,ii-1)+slkm(jj,ll,ii-1);
                        val1 = [val1;tmp];
                    end                    
                end
            end
            %}
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
%                         tmp = 2*pi*auxm(jj,ll,ii-1)*sqrt(P(1,1,ll,ii)*P(2,2,ll,ii)...
%                             -P(1,2,ll,ii)*P(2,1,ll,ii))-this.wt(ll)*exp(-auxgau(jj,ll,ii-1));
                        tmp = 4*pi^2*auxm(jj,ll,ii-1)^2*(P(1,1,ll,ii)*P(2,2,ll,ii)...
                            -P(1,2,ll,ii)*P(2,1,ll,ii))-this.wt(ll)^2*exp(-2*auxgau(jj,ll,ii-1));     
                        val1 = [val1;tmp];
                    end
                end
            end
            
            % auxgau = (x_jj-x_ll)'*Pinv(ll)*(x_jj-x_ll)/2
            for ii = 2:N+1
                for jj = 1:this.gmm_num
                    for ll = 1:this.gmm_num
                        val3 = [val3;auxgau(jj,ll,ii-1)-(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))'*...
                            Pinverse(:,:,ll,ii)*(x(2*jj-1:2*jj,ii)-x(2*ll-1:2*ll,ii))/2];
                    end
                end
            end
           
            %% exp(-auxt(jj)/wt(jj))-sum(auxm(jj,:)) <= 0
%             for ii = 1:N
%                 for jj = 1:this.gmm_num
%                     % log form
%                     val2 = [val2;exp(-auxt(jj,ii)/this.wt(jj))-sum(auxm(jj,:,ii))];
%                     % exp form
%                     val2 = [val2;(-auxt(jj,ii)/this.wt(jj))-log(sum(auxm(jj,:,ii)))+slkt(jj,ii)];
%                 end
%             end
            
            for ii = 1:N
                for jj = 1:this.gmm_num
                    % exp form
                    val2 = [val2;exp(-auxt(jj,ii)/this.wt(jj))-sum(auxm(jj,:,ii))+slkt(jj,ii)];
                    % log form
%                     val2 = [val2;(-auxt(jj,ii)/this.wt(jj))-log(sum(auxm(jj,:,ii)))+slkt(jj,ii)];
                end
            end
            
            h = [val1;val2;val3];
        end
        
        % linear equality constraint
        function val = getLinEqConstr(this,s,snum,tar)
            % linear equality constraints
            z = this.convState(s,snum,'z');
            x = this.convState(s,snum,'x');
            x_pred = this.convState(s,snum,'x_pred');
            P = this.convState(s,snum,'P'); 
            P_pred = this.convState(s,snum,'P_pred'); 
            Pinverse = this.convState(s,snum,'Pinverse');
            
            f = tar.f;
            del_f = tar.del_f;
            Q = tar.Q;
            
            N = this.mpc_hor;
            val = [];
            val1 = [];
            val2 = [];
            val3 = [];
            val4 = [];
            
            % initial condition
            % z(:,1) == this.state;
            % x(:,1) == this.est_pos(:);
            % for jj = 1:this.gmm_num
              %  triu(P(:,:,jj,1)) == triu(this.P{jj});
            % end
            val = [val;z(:,1)-this.state];
            val = [val;x(:,1)-this.est_pos(:)];
            for jj = 1:this.gmm_num
               tmp = triu(P(:,:,jj,1))-triu(this.P{jj}); %%% note: triu here returns a 2*2 matrix, with lower-left item being 0
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
                    % x_k+1|k+1 = x_k+1|k
                    val3 = [val3;x(2*jj-1:2*jj,ii+1)-x_pred(2*jj-1:2*jj,ii)];                    
                end
            end
            
            % symmetric constraint of P and Pinv                                        
            for ii = 1:N+1
                for jj = 1:this.gmm_num
                    val4 = [val4;P(1,2,jj,ii)-P(2,1,jj,ii);Pinverse(1,2,jj,ii)-Pinverse(2,1,jj,ii)];
                end
            end
            
            val = [val;val1;val2;val3;val4];
            %}
        end
        
        % linear inequality constraint
        function glin = getLinIneqConstr(this,s,snum,fld)
            % linear inequality constraints
            z = this.convState(s,snum,'z');
            u = this.convState(s,snum,'u');
            P = this.convState(s,snum,'P'); 
            P_pred = this.convState(s,snum,'P_pred');
            auxm = this.convState(s,snum,'auxm');
            auxt = this.convState(s,snum,'auxt');
            slkt = this.convState(s,snum,'slkt');
            slkm = this.convState(s,snum,'slkm');
            
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
            
            % a weaker constraint for P, Ppred to be psd: P(1,1), P(2,2),
            % P_pred(1,1),P_pred(2,2) >= smnum
            %{
            smnum = 10^-6; % a small number
            for ii = 1:N
                for jj = 1:this.gmm_num
                    glin = [glin;smnum-[P(1,1,jj,ii+1);P(2,2,jj,ii+1);...
                        P_pred(1,1,jj,ii);P_pred(2,2,jj,ii)]];
                end
            end
            %}
            
            % auxiliary variable t,m are nonnegative
            glin = [glin;-auxt(:);-auxm(:)];
            
            %(9/27) try bounding t from above
%             glin = [glin;auxt(:)-10;auxm(:)-1];
%             glin = [glin;slkt(:)-1;slkm(:)-1];
%             glin = [glin;-slkt(:)-1;-slkm(:)-1];
        end
        
        % other convex constraints that are not linear equ/inequ
        % constraints.
        % not finished, not working, and thus not in use.
        % the sdp constraint is kinda hard to formualte in this function to
        % make it usable for the cvx syntax
        function cvxcstr = getCvxConstr(this,s,snum)           
            % psd constraints of P
            P = this.convState(s,snum,'P');
            P_pred = this.convState(s,snum,'P_pred');
            
            cvxcstr = [];
            
            N = this.mpc_hor;
            for ii = 1:N
                for jj = 1:this.gmm_num
                    cvxcstr = [cvxcstr;P_pred(:,:,jj,ii) == semidefinite(2)];
                    cvxcstr = [cvxcstr;P(:,:,jj,ii) == semidefinite(2)];
                end
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
                    tmp = (K(:,:,jj,ii)*(C(:,:,jj,ii)*P_pred(:,:,jj,ii)*C(:,:,jj,ii)'+R)-...
                        P_pred(:,:,jj,ii)*C(:,:,jj,ii)'); % used a scaling factor 
                    val1 = [val1;tmp(:)];
                                        
                    % covariance
                    % P_k+1|k+1=P_k+1|k-gamma*K*C*P_k+1|k
                    tmp_sum = zeros(2,2);
                    T = K(:,:,jj,ii)*C(:,:,jj,ii);
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
               case 'Pinverse'
                   t = s(snum.idx(8,1):snum.idx(8,2));
                   t = reshape(t,2,2,this.gmm_num,N+1);
               case 'auxt'
                   t = s(snum.idx(9,1):snum.idx(9,2));
                   t = reshape(t,this.gmm_num,N);
               case 'auxm'
                   t = s(snum.idx(10,1):snum.idx(10,2));
                   t = reshape(t,this.gmm_num,this.gmm_num,N);
               case 'gam_var'
                   t = s(snum.idx(11,1):snum.idx(11,2));
                   t = reshape(t,this.gmm_num,N);
               case 'slkm'
                   t = s(snum.idx(12,1):snum.idx(12,2));
                   t = reshape(t,this.gmm_num,this.gmm_num,N);
               case 'slkt'
                   t = s(snum.idx(13,1):snum.idx(13,2));
                   t = reshape(t,this.gmm_num,N);
               case 'auxgau'
                   t = s(snum.idx(14,1):snum.idx(14,2));
                   t = reshape(t,this.gmm_num,this.gmm_num,N);
           end
        end
        
        function [s,snum] = setState(this,z,u,x,xpred,P,P_pred,K,Pinverse,auxt,auxm,gam_var,slkm,slkt,auxgau)
            % s is the collection of all decision variables
            % s = [z(:),u(:),x(:),xpred(:),P(:),P_pred(:),K(:)]               
            % when adding new variables, I need to update following
            % functions: setState, convState, genInitState, updOptVar,
            % cvxPlanner_ipopt
            % and the file: labelResult3.m
            s = [z(:);u(:);x(:);xpred(:);P(:);P_pred(:);K(:);Pinverse(:);...
                auxt(:);auxm(:);gam_var(:);slkm(:);slkt(:);auxgau(:)];
            
            % size of each decision variables
            snum = struct();
            snum.z = length(z(:)); %4*(N+1);
            snum.u = length(u(:)); %2*N;
            snum.x = length(x(:)); %2*this.gmm_num*(N+1);
            snum.xpred = length(xpred(:)); %2*this.gmm_num*N;
            snum.P = length(P(:)); %2*2*this.gmm_num*(N+1);
            snum.P_pred = length(P_pred(:)); %2*2*this.gmm_num*N;
            snum.K = length(K(:)); %2*2*this.gmm_num*N;
            snum.Pinverse = length(Pinverse(:)); %2*2*this.gmm_num*N;
            snum.auxt = length(auxt(:)); %this.gmm_num*N
            snum.auxm = length(auxm(:)); %this.gmm_num*this.gmm_num*N
            snum.gamvar = length(gam_var(:)); %this.gmm_num*N
            snum.slkm = length(slkm(:)); %this.gmm_num*this.gmm_num*N
            snum.slkt = length(slkt(:)); %this.gmm_num*N
            snum.auxgau = length(auxgau(:)); %this.gmm_num*this.gmm_num*N
            
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
            snum.idx(8,:) = [snum.idx(7,2)+1,snum.idx(7,2)+snum.Pinverse];
            snum.idx(9,:) = [snum.idx(8,2)+1,snum.idx(8,2)+snum.auxt];
            snum.idx(10,:) = [snum.idx(9,2)+1,snum.idx(9,2)+snum.auxm];
            snum.idx(11,:) = [snum.idx(10,2)+1,snum.idx(10,2)+snum.gamvar];
            snum.idx(12,:) = [snum.idx(11,2)+1,snum.idx(11,2)+snum.slkm];
            snum.idx(13,:) = [snum.idx(12,2)+1,snum.idx(12,2)+snum.slkt];
            snum.idx(14,:) = [snum.idx(13,2)+1,snum.idx(13,2)+snum.auxgau];
            snum.total = snum.idx(14,2);
        end
        
        %% numerical gradient, Jacobian, and hessian
        function [grad,hess] = getNumGradHess(this,f,s,full_hessian)
            % modified from Pieter Abbeel's CS287
%             z = [x(:);P(:)];
%             f = @this.cmpObj;
            y = f(s);
            assert(length(y)==1);           
            
%             fprintf('          [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
%             fprintf('          time of gradient computation.\n');

%             tic
%             grad = gradest(f,s);
%             toc
            
            %
%             tic
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
%                         hess(i,i) = (yhi + ylo - 2*y)/(eps.^2 / 4);
                        grad(i) = (yhi - ylo) / eps;
                    end
                else
                    grad = this.numerical_jac(f,s);
                    hess = this.numerical_jac(@(z) this.numerical_jac(f,z), s);
                    hess = (hess + hess')/2;
                end
%                 mineigv = min(eig(hess));
%                 hess = hess+(mineigv+0.1)*eye(length(s));
            end
%             toc
            %}
        end
        
        % jacobian
        function jaco = getNumJac(this,f,s)
            
%             fprintf('          [inner sqp loop] Robot.m, line %d\n', MFileLineNr())
%             fprintf('          time of Jacobian computation.\n');

%             tic
%             jaco = jacobianest(f,s);
%             toc

            % modified from Pieter Abbeel's CS287            
            %
%             tic
            y = f(s);
            
            jaco = zeros(length(y), length(s));
            
            eps = 1e-5;
            sp = s;
            
            for i=1:length(s)
                sp(i) = s(i) + eps/2;
                yhi = f(sp);
                sp(i) = s(i) - eps/2;
                ylo = f(sp);
                sp(i) = s(i);
                jaco(:,i) = (yhi - ylo) / eps;
            end  
%             toc
            %}
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
        
        function gam_exact = gam_eps(this,z,theta,x0,alp1,alp2,alp3)
            % gamma = (eps+exp^x)/(1+exp^x)
            eps1 = 10^-3;
            [~,pow1] = this.gam_den1(z,x0,alp1);
            [~,pow2] = this.gam_den2(z,x0,theta,alp2);
            [~,pow3] = this.gam_den3(z,x0,theta,alp3);
            gam_exact = ((eps1+exp(-pow1))*(eps1+exp(-pow2))*(eps1+exp(-pow3)))/...
                ((1+exp(-pow1))*(1+exp(-pow2))*(1+exp(-pow3)));
        end
                
        % denominator of gamma
        function gam_d = gam_den(this,z,theta,x0,alp1,alp2,alp3)
            [gam_d1,~] = this.gam_den1(z,x0,alp1);
            [gam_d2,~] = this.gam_den2(z,x0,theta,alp2);
            [gam_d3,~] = this.gam_den3(z,x0,theta,alp3);
            gam_d = gam_d1*gam_d2*gam_d3;
%             gam_d = this.gam_den1(z,x0,alp1)*this.gam_den2(z,x0,theta,alp2)*...
%                 this.gam_den3(z,x0,theta,alp3);
        end
        
        function [gam_den,gam_pow] = gam_den1(this,z,x0,alp)
            tmp = alp*(sum((x0-z).^2)-this.r^2);
            % yalmip cannot accept if statement in constraint
            gam_pow = tmp;
            gam_den = 1+exp(tmp);            
        end
        
        function [gam_den,gam_pow] = gam_den2(this,z,x0,theta,alp)
            tmp = alp*[sin(theta-this.theta0),-cos(theta-this.theta0)]*(x0-z);
            gam_pow = tmp;
            gam_den = 1+exp(tmp);            
        end
        
        function [gam_den,gam_pow] = gam_den3(this,z,x0,theta,alp)
            tmp = alp*[-sin(theta+this.theta0),cos(theta+this.theta0)]*(x0-z);
            gam_pow = tmp;
            gam_den = 1+exp(tmp);            
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
            % Plotting FOV along the planned path
            
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