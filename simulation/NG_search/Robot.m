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
        gam; % function handle for gamma
        gam_den; % function handle for the denominator of gamma
        gam_aprx; % function handle for the linearized gamma
        p_aprx; % function handle for the element of the linearized covriance matrix
        alp1; % parameters in gam
        alp2; % parameters in gam
        alp3; % parameters in gam
        
        % observation
        y; % observation measurement
        
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
            this.gam = inPara.gam; % function handle for gamma
            this.gam_den = inPara.gam_den; % function handle for the denominator of gamma
            this.gam_aprx = inPara.gam_aprx; % function handle for the linearized gamma
            this.p_aprx = inPara.p_aprx; % function handle for the element of the linearized covriance matrix
            this.alp1 = inPara.alp1; % parameters in gam
            this.alp2 = inPara.alp2;
            this.alp3 = inPara.alp3;
            
            % filtering
            this.sensor_type = inPara.sensor_type;
            % xKF
            this.est_pos = inPara.est_pos;
            this.P = inPara.P;
            this.est_pos_hist = [];
            this.P_hist = [];
            % gmm
            this.gmm_num = inPara.gmm_num;
            this.wt = inPara.wt;
            % pf
            this.max_gmm_num = inPara.max_gmm_num;
            this.particles = inPara.particles;
            
            this.mpc_hor = inPara.mpc_hor;
            this.dt = inPara.dt;
            this.optu = [];
            this.gam = inPara.gam;
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
            %             flag = (norm(tar_pos-this.state(1:2)) <= this.r);
        end
        
        %% measurement generation
        % generate a random measurement
        function y = sensorGen(this,fld)
            tar_pos = fld.target.pos;
            % range-bearing sensor
            if strcmp(this.sensor_type,'rb')
                if this.inFOV(tar_pos)
                    y = this.h(tar_pos,this.state(1:2))+(mvnrnd([0;0],this.R))';
                else
                    y = [-100;-100];
                end
            elseif strcmp(this.sensor_type,'ran')
                if this.inFOV(tar_pos)
                    y = norm(tar_pos-this.state(1:2))+normrnd(0,this.R);
                else
                    y = -100;
                end
            end
        end
        
        %% filtering
        function this = KF(this,fld)
            y = this.y;
            
            % target
            tar = fld.target;
            A = tar.A;
            B = tar.B;
            Q = tar.Q;
            
            % current estimation
            x = this.est_pos;
            P = this.P;
            
            % sensor
            C = this.C;
            R = this.R;
            
            % prediction
            x_pred = A*x+B;
            P_pred = A*P*A'+Q;
            
            % update
            if sum(y-[-100;-100]) ~= 0
                % if an observation is obtained
                K = P_pred*C'/(C*P_pred*C'+R);
                x_next = x_pred+K*(y-C*x_pred);
                P_next = P_pred-K*C*P_pred;
            else
                x_next = x_pred;
                P_next = P_pred;
            end
            
            this.est_pos = x_next;
            this.P = P_next;
            this.est_pos_hist = [this.est_pos_hist,x_next];
            this.P_hist = [this.P_hist,P_next];
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
                        w(ii) = 0.01;
                    else
                        w(ii) = 0.99;
                    end
                else
                    if this.inFOV(pred_par(:,ii))
                        w(ii) = mvnpdf(y,this.h(pred_par(:,ii),this.state(1:2)),this.R);
                    else
                        w(ii) = 0.01;
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
            display ('gmm fitting takes time as:');
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
        
        % try re-writing the problem using cvx solver
        function [optz,optu] = cvxPlanner(this,fld,optz,optu) % cvxPlanner(this,fld,init_sol)
            % use the multi-layer approach similar to Sachin's work. Fix
            % the parameter for the sensor, solve path planning. Then
            % refine the parameter until close to reality. In each
            % iteration, a convex program is solved. The initial solution
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
            gam = this.gam;
            gam_aprx = this.gam_aprx;
           	p_aprx = this.p_aprx;
            
            % the parameter for the sensing boundary approximation
            alp = 1;
            alp_inc = 2; % increament paramter for alpha
            
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
            P_pred_ref = init_sol.P_pred;
            
%             zref = init_sol.zref;
%             uref = init_sol.uref;
%             xref = init_sol.xref;
%             Kref = init_sol.Kref;
%             P_pred_ref = init_sol.P_pred_ref;
            
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
                
                % obj
                minimize sum(t(:))
                               
                % constraints
                % epigraph for obj
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
                        end
                    end
                end
                
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
                    
                    % use the weighted mean as the MAP of target position
                    
%                     tmp_mean = reshape(xref(:,ii+1),2,this.gmm_num)*this.wt;
                    
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
                        
%                         triu(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii)) == -triu(tmp);
%                         triu(P(:,:,jj,ii+1)-P_pred(:,:,jj,ii))*gamma_den...
%                             == -gamma_num*triu(Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C*P_pred(:,:,jj,ii));
                    end
                end
                this.w_lb <= u(1,:) <= this.w_ub;
                this.a_lb <= u(2,:) <= this.a_ub;
                this.v_lb <= z(4,:) <= this.v_ub;
                
                cvx_end
                
                % terminating condition: the actual in/out FOV is
                % consistent with that of planning
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
                
                dif = max(max(abs(is_in_fov-gamma_aprx)));
                if dif < 0.05
                    break
                end               
                
                % assign value for next iteration
                zref = z;
                uref = u;
                xref = x;                
                P_pred_ref = P_pred;
                
                % coompute Kref using Ricatti equation 
                for ii = 1:N
                    for jj = 1:this.gmm_num
                        C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
                        Kref(2*jj-1:2*jj,2*ii-1:2*ii) = P_pred(1,2,jj,ii)*C'/(C*P_pred(1,2,jj,ii)*C'+R);                        
                    end
                end
                
                
                alp1 = alp1*alp_inc;
                alp2 = alp2*alp_inc;
                alp3 = alp3*alp_inc;
                
                display ('Robot.m line 522')
                display (alp1)
                display (alp2)
                display (alp3)
            end
            
            optz = zref;
            optu = uref;
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
            % generate initial states for ngPlanner
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
                    gamma = this.inFOV(x_olp(2*jj-1:2*jj,ii));
                    P(:,:,jj,ii+1) = P_pred(:,:,jj,ii)-gamma*K(2*jj-1:2*jj,2*ii-1:2*ii)*P_pred(:,:,jj,ii)*C;
                end
            end
            init_state.P = P;
            init_state.P_pred = P_pred;
            init_state.K = K;            
        end
               
        %% robot state updating
        function this = updState(this,u)
            st = this.state;
            this.optu = [this.optu,u(:,1)];
            dt = this.dt;
            this.state = st+[st(4)*cos(st(3));st(4)*sin(st(3));u(:,1)]*dt;
            this.traj = [this.traj,this.state];
            % range-bearing sensor
%             this.h = @(x) x-this.state(1:2);
        end
    end
end