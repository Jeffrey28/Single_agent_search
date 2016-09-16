classdef Robot
    properties
        % motion specs
        traj;
        state; % [x;y;theta;v]
        v_lb;
        v_ub;
        w_lb;
        w_ub;
        
        % sensor specs
        theta0; % sensor range in angle
        r; % sensor range
        C; % C matrix of observation
        R; % covariance for sensor model
        
        % observation
        y; % observation measurement
        
        % filtering
        est_pos; % estimated target position
        P; % estimation covariance
        
        % path planning
        mpc_hor;
        dt;
        optu;
        gam; % coefficient for sigmoid function
        
        % performance metrics
        ml_pos;
        ent_pos;
    end
    
    methods
        function this = Robot(inPara)
            this.state = inPara.state;
            this.traj = inPara.state(1:2);
            this.v_lb = inPara.v_lb;
            this.v_ub = inPara.v_ub;
            this.w_lb = inPara.w_lb;
            this.w_ub = inPara.w_ub;
            this.R = inPara.R;
            this.C = inPara.C;
            this.theta0 = inPara.theta0;
            this.r = inPara.range;
            this.est_pos = inPara.init_pos;
            this.P = inPara.init_P;
            
            this.mpc_hor = inPara.mpc_hor;
            this.dt = inPara.dt;
            this.optu = [];
            this.gam = inPara.gam;
        end
        
        % approximate straightline edge of sensor FOV based on current
        % robot state
        function [a,b] = FOV(this,st)
            theta = st(3);
            x0 = st(1);
            y0 = st(2);
            alp1 = theta - this.theta0;
            alp2 = theta + this.theta0;
            a = [sin(alp1),-cos(alp1);-sin(alp2),cos(alp2)]; % [a1;a2]
            b = [-x0*sin(alp1)+y0*cos(alp1);x0*sin(alp2)-y0*cos(alp2)];%[b1;b2];
        end
        
        % determine if the target is in sensor FOV
        function flag = inFOV(this,tar_pos)
            [a,b] = this.FOV(this.state);
            flag = (a(1,:)*tar_pos-b(1) <= 0) && (a(2,:)*tar_pos-b(2) <= 0)...
                && (norm(tar_pos-this.state(1:2)) <= this.r);
        end
        
        % generate a random measurement
        function y = sensorGen(this,fld)
            tar_pos = fld.target.pos;
            
            if this.inFOV(tar_pos)
                y = tar_pos+mvnrnd([0;0],this.R);
            else
                y = [-100;-100];
            end
        end
        
        function this = KF(this,fld)
            y = this.y;
            
            % target
            tar = fld.target;
            A = tar.A;
            Q = tar.Q;
            
            % current estimation
            x = this.est_pos;
            P = this.P;
            
            % sensor
            C = this.C;
            R = this.R;
            
            % prediction
            x_pred = A*x;
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
        end
        
        function [optz,optu] = Planner(this,fld)
            N = this.mpc_hor;
            var_dt = this.dt;
            var_C = this.C;
            var_R = this.R;
            st = this.state;
            x_cur = this.est_pos;
            gam = this.gam;
            
            tar = fld.target;
            A = tar.A;
            Q = tar.Q;
            
            % set up simulation
            % robot state and control
            z = sdpvar(4,N+1,'full');
            u = sdpvar(2,N,'full');
            % estimation
            x = sdpvar(2,N+1,'full');
%             P = sdpvar(2,2,N+1);
%             var_P = sdpvar(2,2,N+1,'full'); % symmetric matrix in first two dim
            var_P = sdpvar(2,2*(N+1),'full'); % symmetric matrix in first two dim
            % auxiliary variable
            tmp_M = sdpvar(2,2,'full');
            K = sdpvar(2,2,'full');
            phi = sdpvar(2,2,'full');
            
            % debug purpose
            x_pred = sdpvar(2,N,'full');
%             P_pred = sdpvar(2,2,N,'full');
            P_pred = sdpvar(2,2*N,'full');
            
            % obj
            obj = 0;%P(1,1,N+1)+P(2,2,N+1); % trace of last covariance
            
            % constraints
            constr = [z(:,1) == this.state];
            constr = [constr,x(:,1) == this.est_pos];
            constr = [constr,var_P(:,1:2) == this.P];%[1 0;0 1]];
%             constr = [constr,var_P(:,:,1) == [1 0;0 1]];%this.P];
            
            for ii = 1:N
                obj = obj+var_P(1,2*ii+1)+var_P(2,2*ii+2);%var_P(1,1,ii+1)^2+var_P(2,2,ii+1)^2+sum(sum(phi.^2));
                
%                 constr = [constr,P(:,:,ii+1)>=0];
                
                % robot state
                constr = [constr,z(:,ii+1) == z(:,ii)+...
                    [z(4,ii)*cos(z(3,ii));z(4,ii)*sin(z(3,ii));...
                    u(:,ii)]*var_dt];
                
                constr = [constr,[fld.fld_cor(1);fld.fld_cor(3)]<=z(1:2,ii+1)<=...
                    [fld.fld_cor(2);fld.fld_cor(4)]];
                
                % KF update
%                 alp1 = z(3,ii+1) - this.theta0;
%                 alp2 = z(3,ii+1) + this.theta0;
                alp1 = st(3) - this.theta0;
                alp2 = st(3) + this.theta0;
                a = [sin(alp1),-cos(alp1);-sin(alp2),cos(alp2)]; % [a1;a2]
                b = [-z(1,ii+1)*sin(alp1)+z(2,ii+1)*cos(alp1);z(1,ii+1)*sin(alp2)-z(2,ii+1)*cos(alp2)];%[b1;b2];
                
%                 alp1 = st(3) - this.theta0;
%                 alp2 = st(3) + this.theta0;
%                 a = [sin(alp1),-cos(alp1);-sin(alp2),cos(alp2)]; % [a1;a2]
%                 b = [-st(1)*sin(alp1)+st(2)*cos(alp1);st(1)*sin(alp2)-st(2)*cos(alp2)];%[b1;b2];
                
%                 delta = 1;%/((1+exp(gam*(a(1,:)*x(:,ii)-b(1))))*(1+exp(gam*(a(2,:)*x(:,ii)-b(2))))*...
%                     (1+exp(gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2))));
%                 delta_rcp = (1+exp(gam*(a(1,:)*x(:,ii)-b(1))))*(1+exp(gam*(a(2,:)*x(:,ii)-b(2))))*...
%                     (1+exp(gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2)));
%                 D = ...%(1+exp(gam*(a(1,:)*x(:,ii)-b(1))))*(1+exp(gam*(a(2,:)*x(:,ii)-b(2))))*...
%                     (1+exp(gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2)));
%                 dist_prod = (gam*(a(1,:)*x(:,ii)-b(1)))^2;%(gam*(a(1,:)*x(:,ii)-b(1))*gam*(a(2,:)*x(:,ii)-b(2))*...
% %                     gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2))^2;
%                 dist_prod2 = (1+(gam*(a(1,:)*x(:,ii)-b(1)))^2);%*(1+(gam*(a(2,:)*x(:,ii)-b(2)))^2)*...
% %                     (1+(gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2))^2);
                dist_prod = 1;% (sqrtm(sum((x(:,ii)-this.state(1:2))).^2)-this.r)^2;%(gam*(a(1,:)*x_cur-b(1)))^2*(gam*(a(2,:)*x(:,ii)-b(2)))^2;%(gam*(a(1,:)*x(:,ii)-b(1))*gam*(a(2,:)*x(:,ii)-b(2))*...
%                     gam*(sum((x(:,ii)-this.state(1:2)).^2)-this.r^2))^2;
                dist_prod2 = ...%(1+(gam*(a(1,:)*x_cur-b(1)))^2)*(1+(gam*(a(2,:)*x(:,ii)-b(2)))^2);%*...
                    1+sum((x(:,ii)-this.state(1:2)).^2);%(1+(sqrtm(sum((x(:,ii)-this.state(1:2)).^2))-this.r)^2);
                

                % prediction
                x_pred(:,ii) = A*x(:,ii);
%                 P_pred(:,:,ii) = A*var_P(:,:,ii)*A'+Q;
                P_pred(:,2*ii-1:2*ii) = A*var_P(:,2*ii-1:2*ii)*A'+Q;
                
                % update
% %                 K = P_pred*C'*delta*tmp_M*delta;
% %                 constr = [constr,K == 1];
% %                 constr = [constr,K*(var_C*P_pred(:,:,ii)*var_C'+var_R/delta^2)==P_pred(:,:,ii)*var_C'];
%                 constr = [constr,K*(var_C*P_pred(:,2*ii-1:2*ii)*var_C'+var_R/delta^2)==P_pred(:,2*ii-1:2*ii)*var_C'];
% %                 constr = [constr,K*(C*P_pred(:,:,ii)*C'+R*D^2)==P_pred(:,:,ii)*C'];
%                 constr = [constr,x(:,ii+1) == x_pred(:,ii)+K*(var_C*x(:,ii)-var_C*x_pred(:,ii))];
% %                 constr = [constr,x(:,ii+1) == x_pred+K*(C*x_cur-C*x_pred)];
% %                 constr = [constr,var_P(:,:,ii+1) == P_pred(:,:,ii)-K*var_C*P_pred(:,:,ii)+phi];
%                 constr = [constr,var_P(:,2*ii+1:2*ii+2) >= P_pred(:,2*ii-1:2*ii)-K*var_C*P_pred(:,2*ii-1:2*ii)];%+phi];
% %                 constr = [constr,(delta*C*P_pred(:,:,ii)*C'*delta+R)*tmp_M == eye(2)];
                
%                 % test version of simplified constraint
%                 T = var_C*P_pred(:,2*ii-1:2*ii)*var_C'+var_R*delta_rcp^2;
%                 a = T(1,1);
%                 b = T(1,2);
%                 c = T(2,1);
%                 d = T(2,2);
%                 t = a*d-b*c;
%                 T2 = [d -b; -c a];
%                 constr = [constr,(x(:,ii+1)-x_pred(:,ii))*t == P_pred(:,2*ii-1:2*ii)*var_C'*T2*(var_C*x(:,ii)-var_C*x_pred(:,ii))];
%                 constr = [constr,(var_P(:,2*ii+1:2*ii+2)-P_pred(:,2*ii-1:2*ii))*t...
%                     == -P_pred(:,2*ii-1:2*ii)*var_C'*T2*var_C*P_pred(:,2*ii-1:2*ii)];%+phi];
                
                % use x/sqrt(1+x^2) as sigmoid function and Sachin's
                % formulation
                T = var_C*P_pred(:,2*ii-1:2*ii)*var_C'*dist_prod+var_R*dist_prod2;
                a = T(1,1);
                b = T(1,2);
                c = T(2,1);
                d = T(2,2);
                t = a*d-b*c;
                T2 = [d -b; -c a];
                constr = [constr,(x(:,ii+1)-x_pred(:,ii))*t == dist_prod*P_pred(:,2*ii-1:2*ii)*var_C'*T2*(var_C*x(:,ii)-var_C*x_pred(:,ii))];
                constr = [constr,(var_P(:,2*ii+1:2*ii+2)-P_pred(:,2*ii-1:2*ii))*t...
                    == -dist_prod*P_pred(:,2*ii-1:2*ii)*var_C'*T2*var_C*P_pred(:,2*ii-1:2*ii)];%+phi];
                
                % use Schenato's formulation
%                 T = var_C*P_pred(:,2*ii-1:2*ii)*var_C'*dist_prod+var_R;
%                 a = T(1,1);
%                 b = T(1,2);
%                 c = T(2,1);
%                 d = T(2,2);
%                 t = a*d-b*c;
%                 T2 = [d -b; -c a];
%                 constr = [constr,(x(:,ii+1)-x_pred(:,ii))*dist_prod2*t == dist_prod*P_pred(:,2*ii-1:2*ii)*var_C'*T2*(var_C*x(:,ii)-var_C*x_pred(:,ii))];
%                 constr = [constr,(var_P(:,2*ii+1:2*ii+2)-P_pred(:,2*ii-1:2*ii))*dist_prod2*t...
%                     == -dist_prod*P_pred(:,2*ii-1:2*ii)*var_C'*T2*var_C*P_pred(:,2*ii-1:2*ii)];%+phi];
            end
            constr = [constr, this.w_lb <= u(1,:) <= this.w_ub, this.v_lb <= u(2,:) <= this.v_ub];
            
            opt = sdpsettings('solver','ipopt','verbose',3,'debug',1,'showprogress',1);
            
            sol = optimize(constr,obj,opt);
            optz = value(z);
            optu = value(u);
        end
        
        function this = updState(this,u)
            st = this.state;
            this.optu = [this.optu,u(:,1)];
            dt = this.dt;
            this.state = st+[st(4)*cos(st(3));st(4)*sin(st(3));u(:,1)]*dt;
        end
        
        function this = computeMetrics(this,fld,id)
            % Computing Performance Metrics
            count = this.step_cnt;
            
            %% ML error
            % DBF
            if strcmp(id,'dbf')
                [tmp_x1,tmp_y1] = find(this.dbf_map == max(this.dbf_map(:)));
                if length(tmp_x1) > 1
                    tmp_idx = randi(length(tmp_x1),1,1);
                else
                    tmp_idx = 1;
                end
                this.ml_pos_dbf(:,count) = [tmp_x1(tmp_idx);tmp_y1(tmp_idx)];
                this.ml_err_dbf(count) = norm(this.ml_pos_dbf(:,count)-fld.target.pos);
            end
            
            % concensus
            if strcmp(id,'cons')
                [tmp_x2,tmp_y2] = find(this.cons_map == max(this.cons_map(:)));
                if length(tmp_x2) > 1
                    tmp_idx2 = randi(length(tmp_x2),1,1);
                else
                    tmp_idx2 = 1;
                end
                this.ml_pos_cons(:,count) = [tmp_x2(tmp_idx2);tmp_y2(tmp_idx2)];
                this.ml_err_cons(count) = norm(this.ml_pos_cons(:,count)-fld.target.pos);
            end
            
            % centralized
            if strcmp(id,'cent')
                [tmp_x3,tmp_y3] = find(this.cent_map == max(this.cent_map(:)));
                if length(tmp_x3) > 1
                    tmp_idx3 = randi(length(tmp_x3),1,1);
                else
                    tmp_idx3 = 1;
                end
                this.ml_pos_cent(:,count) = [tmp_x3(tmp_idx3);tmp_y3(tmp_idx3)];
                this.ml_err_cent(count) = norm(this.ml_pos_cent(:,count)-fld.target.pos);
            end
            
            %% Covariance of posterior pdf
            [ptx,pty] = meshgrid(1:fld.fld_size(1),1:fld.fld_size(2));
            pt = [ptx(:),pty(:)];
            
            % DBF
            if strcmp(id,'dbf')
                tmp_map1 = this.dbf_map;
                % this avoids the error when some grid has zeros probability
                tmp_map1(tmp_map1 <= realmin) = realmin;
                
                % compute covariance of distribution
                dif1 = pt' - [(1+fld.target.pos(1))/2;(1+fld.target.pos(2))/2]*ones(1,size(pt',2));
                cov_p1 = zeros(2,2);
                for jj = 1:size(pt',2)
                    cov_p1 = cov_p1 + dif1(:,jj)*dif1(:,jj)'*tmp_map1(pt(jj,1),pt(jj,2));
                end
                this.pdf_cov_dbf{count} = cov_p1;
                this.pdf_norm_dbf(count) = norm(cov_p1,'fro');
            end
            
            if strcmp(id,'cons')
                % concensus
                tmp_map2 = this.cons_map;
                % this avoids the error when some grid has zeros probability
                tmp_map2(tmp_map2 <= realmin) = realmin;
                
                % compute covariance of distribution
                dif2 = pt' - [(1+fld.target.pos(1))/2;(1+fld.target.pos(2))/2]*ones(1,size(pt',2));
                cov_p2 = zeros(2,2);
                for jj = 1:size(pt',2)
                    cov_p2 = cov_p2 + dif2(:,jj)*dif2(:,jj)'*tmp_map2(pt(jj,1),pt(jj,2));
                end
                this.pdf_cov_cons{count} = cov_p2;
                this.pdf_norm_cons(count) = norm(cov_p2,'fro');
            end
            
            % centralized
            if strcmp(id,'cent')
                tmp_map3 = this.cent_map;
                % this avoids the error when some grid has zeros probability
                tmp_map3(tmp_map3 <= realmin) = realmin;
                
                % compute covariance of distribution
                dif3 = pt' - [(1+fld.target.pos(1))/2;(1+fld.target.pos(2))/2]*ones(1,size(pt',2));
                cov_p3 = zeros(2,2);
                for jj = 1:size(pt',2)
                    cov_p3 = cov_p3 + dif3(:,jj)*dif3(:,jj)'*tmp_map3(pt(jj,1),pt(jj,2));
                end
                this.pdf_cov_cent{count} = cov_p3;
                this.pdf_norm_cent(count) = norm(cov_p3,'fro');
            end
            
            %% Entropy of posterior pdf
            % DBF
            if strcmp(id,'dbf')
                tmp_map1 = this.dbf_map;
                % this avoids the error when some grid has zeros probability
                tmp_map1(tmp_map1 <= realmin) = realmin;
                dis_entropy = -(tmp_map1).*log2(tmp_map1); % get the p*log(p) for all grid points
                this.ent_dbf(count) = sum(sum(dis_entropy));
            end
            
            % concensus
            if strcmp(id,'cons')
                tmp_map2 = this.cons_map;
                % this avoids the error when some grid has zeros probability
                tmp_map2(tmp_map2 <= realmin) = realmin;
                dis_entropy = -(tmp_map2).*log2(tmp_map2); % get the p*log(p) for all grid points
                this.ent_cons(count) = sum(sum(dis_entropy));
            end
            
            % centralized
            if strcmp(id,'cent')
                tmp_map3 = this.cent_map;
                % this avoids the error when some grid has zeros probability
                tmp_map3(tmp_map3 <= realmin) = realmin;
                dis_entropy = -(tmp_map3).*log2(tmp_map3); % get the p*log(p) for all grid points
                this.ent_cent(count) = sum(sum(dis_entropy));
            end
        end
    end
end