classdef Robot
    properties        
        % motion state
        traj;
        pos; % robot position
        
        % sensor specs
        sen_cov;
        inv_sen_cov;
        
        % observation
        z; % observation measurement
        
        % filtering
        lkhd_map; % prob matrix for a certain observation. a way to reduce computation time, sacrificing the space complexity
        
        % performance metrics
        ml_pos;
        ent_pos;
    end
    
    methods
        function this = Robot(inPara)
            this.pos = inPara.pos;
            this.traj = inPara.pos;
            this.sen_cov = inPara.sen_cov;
            this.inv_sen_cov = inPara.inv_sen_cov;
%             this.lkhd_map = zeros(inPara.fld_size(1),inPara.fld_size(2));                       
        end
        
        % generate a random measurment and computes the probability
        % likelihood map
        function this = sensorGen(this,fld)
            % generate sensor measurement
            x_r = this.pos;
            x_t = fld.target.pos;
            inv_cov = this.inv_sen_cov;
            offset = this.sen_offset;
            
            tmp_lkhd = exp(-1/2*(x_t+offset-x_r)'*inv_cov*(x_t+offset-x_r));
            tmp_z = (rand(1,1) < tmp_lkhd);
            
            this.z = tmp_z;
            this.k = this.step_cnt;
            
            % generate the likelihood map for all possible target locations
            tmp_lkhd_map = this.sensorProb(fld);
            if tmp_z == 0
                tmp_lkhd_map = 1-tmp_lkhd_map;
            end
            
            this.lkhd_map = tmp_lkhd_map;
        end
        
        function this = KF(this,fld)
            
        end
        
%         % computes probability likelihood map
%         function lkhd_map = sensorProb(this,fld)
%             x_r = this.pos;
%             inv_cov = this.inv_sen_cov;
%             offset = this.sen_offset;
%             
%             xlen = fld.fld_size(1);
%             ylen = fld.fld_size(2);
%             [ptx,pty] = meshgrid(1:xlen,1:ylen);
%             pt = [ptx(:)';pty(:)'];
%             pt = bsxfun(@plus,pt,offset);
%             pt = bsxfun(@minus,pt,x_r);
%             tmp = pt'*inv_cov*pt;
%             tmp_diag = diag(tmp);
%             lkhd_map = exp(-1/2*(reshape(tmp_diag,ylen,xlen))');% note meshgrid first goes along y and then x direction.
%         end
        
%         function this = updOwnMsmt(this,inPara)
%             % (1) self-observation
%             % update robot's own measurement in the communication buffer
%             
%             selection = inPara.selection;
%             if (selection == 1) || (selection == 3)
%                 % static target
%                 % update the buffer with the robot's own observation
%                 this.buffer(this.idx).pos = this.pos;
%                 this.buffer(this.idx).z = this.z;
%                 this.buffer(this.idx).k = this.step_cnt;
%                 this.buffer(this.idx).lkhd_map = this.lkhd_map;
%                 
%                 % assign this probability to rbt_cons and rbt_cent to
%                 % save computation resource
%                 %                 rbt.cons_prob = rbt.buffer(rbt.idx).prob;
%                 %                 rbt.cent_prob = rbt.buffer(rbt.idx).prob;
%                 %                 rbt_cons(i).prob = rbtBuffer{i}.rbt(i).prob;
%                 %                 rbt_cent.prob{i} = rbtBuffer{i}.rbt(i).prob;
%                 
%             elseif (selection == 2) || (selection == 4)
%                 % moving target
%                 this.buffer(this.idx).pos = [this.buffer(this.idx).pos,this.pos];
%                 this.buffer(this.idx).z = [this.buffer(this.idx).z,this.z];
%                 this.buffer(this.idx).k = [this.buffer(this.idx).k,this.k];
%                 this.buffer(this.idx).lkhd_map{this.step_cnt} = this.lkhd_map;
%                 
%                 % assign this probability to rbt_cons and rbt_cent to
%                 % save computation resource
%                 %                 rbt.cons_prob = rbt.buffer(rbt.idx).prob;
%                 %                 rbt.cent_prob = rbt.buffer(rbt.idx).prob;                
%             end
%         end
%         
%         function this = dataExch(this,inPara)
%             % (2) data transmission
%             % exchange communication buffer with neighbors
%             selection = inPara.selection;
%             rbt_nbhd_set = inPara.rbt_nbhd_set;
%             if (selection == 1) || (selection == 3)
%                 % static target
%                 % exchange and update buffer
%                 for t = 1:length(rbt_nbhd_set)
%                     % note: communication only transmit the latest
%                     % observation stored in each neighbor
%                     tmp_rbt = rbt_nbhd_set{t};
%                     for jj = 1:this.num_robot
%                         if (~isempty(tmp_rbt.buffer(jj).k)) && (isempty(this.buffer(jj).k) || (this.buffer(jj).k < tmp_rbt.buffer(jj).k))
%                             this.buffer(jj).pos = tmp_rbt.buffer(jj).pos;
%                             this.buffer(jj).z = tmp_rbt.buffer(jj).z;
%                             this.buffer(jj).k = tmp_rbt.buffer(jj).k;
%                             this.buffer(jj).lkhd_map = tmp_rbt.buffer(jj).lkhd_map;
%                         end
%                     end
%                 end                
%                 
%             elseif (selection == 2) || (selection == 4)
%                 % moving target
%                 for t = 1:length(rbt_nbhd_set)
%                     % note: communication only transmit the latest
%                     % observation stored in each neighbor
%                     tmp_rbt = rbt_nbhd_set{t};
%                     for jj = 1:this.num_robot
%                         if (~isempty(tmp_rbt.buffer(jj).k)) && (isempty(this.buffer(jj).k) || (this.buffer(jj).k(end) < tmp_rbt.buffer(jj).k(end)))
%                             %%% this code only handles the fixed topology
%                             %%% case, i.e. each time a one-step newer
%                             %%% observation is received. for multi-step
%                             %%% newer observations, such code needs
%                             %%% modification.
%                             this.buffer(jj).pos = [this.buffer(jj).pos,tmp_rbt.buffer(jj).pos(:,end)];
%                             this.buffer(jj).z = [this.buffer(jj).z,tmp_rbt.buffer(jj).z(end)];
%                             this.buffer(jj).k = [this.buffer(jj).k,tmp_rbt.buffer(jj).k(end)];
%                             % in real experiment, this prob term should not
%                             % be communicated. In simulation, this is for
%                             % the purpose of accelerating the computation speed.
%                             if isempty(this.buffer(jj).lkhd_map)                                
%                                 this.buffer(jj).lkhd_map = tmp_rbt.buffer(jj).lkhd_map(end);
%                             else
%                                 this.buffer(jj).lkhd_map(end+1) = tmp_rbt.buffer(jj).lkhd_map(end);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         
%         function this = DBF(this,inPara)
%             % filtering
%             selection = inPara.selection;
%             if (selection == 1) || (selection == 3)
%                 % calculate probility of latest z               
%                 for jj=1:this.num_robot % Robot Iteration
%                     if (~isempty(this.buffer(jj).k)) && (~ismember(this.buffer(jj).k,this.buffer(jj).used))
%                         this.dbf_map=this.dbf_map.*this.buffer(jj).lkhd_map;
%                         this.buffer(jj).used = [this.buffer(jj).used,this.buffer(jj).k];
%                     end
%                 end
%                 this.dbf_map=this.dbf_map/sum(sum(this.dbf_map));
%                                 
%             elseif (selection == 2) || (selection == 4)
%                 upd_matrix = this.upd_matrix{1};
%                 %% update by bayes rule
%                 % note: main computation resource are used in calling sensorProb function.
%                 % when using grid map, can consider precomputing
%                 % results and save as a lookup table
%                
%                 talign_flag = 1; % if all agent's observation's time are no less than talign_t+1, then talign_flag = 1, increase talign_t
%                 tmp_t = this.talign_t;
%                 tmp_map = this.talign_map; % time-aligned map
%                 
% %                 display(this.idx)
%                 for t = (this.talign_t+1):this.step_cnt
%                     display(t)
%                     %% one-step prediction step                     
%                     % p(x_k+1) = sum_{x_k} p(x_k+1|x_k)p(x_k)
%                     % note, data in upd_matrix is first along y direction 
%                     % and then along x direction. therefore, when using
%                     % tmp_map(:), care should be taken since tmp_map(:)
%                     % orders data first along x direction, then y
%                     % direction.
%                     tmp_map = tmp_map';
%                     tmp_map2 = upd_matrix*tmp_map(:);
%                     tmp_map = (reshape(tmp_map2,size(tmp_map)))';
%                     
%                     %% updating step                                       
%                     for jj=1:this.num_robot
%                         display(jj)
%                         if (~isempty(this.buffer(jj).k)) && (this.buffer(jj).k(end) >= t)
%                             % note: this update is not valid in real
%                             % experiment since we don't communicate
%                             % probability. This line of code is for
%                             % computation reduction in simulation
%                             tmp_map = tmp_map.*this.buffer(jj).lkhd_map{t};
%                         else
%                             talign_flag = 0;
%                         end
%                     end
%                     
%                     % assign the probability to rbt_cons and rbt_cent to
%                     % save computation resource
%                     %                         rbt_cons(ii).prob = rbtBuffer{ii}.rbt(ii).prob;
%                     %                         rbt_cent.prob{ii} = rbtBuffer{ii}.rbt(ii).prob;
%                     
%                     % after the first loop, the robot's aligned time
%                     % increase by one if the robot's buffer is no later
%                     % than the previous aligned time
%                     if (t == this.talign_t+1) && (talign_flag == 1)
%                         this.talign_map = tmp_map;
%                         this.talign_map = this.talign_map/sum(sum(this.talign_map));
%                         tmp_t = tmp_t+1;
%                         % when we increase the aligned time, we can remove 
%                         % the previous data 
%                         if this.talign_t > 0 
%                             for jj=1:this.num_robot
%                                 this.buffer(jj).pos(:,this.talign_t) = -ones(2,1);
%                                 this.buffer(jj).z(this.talign_t) = -1;
%                                 this.buffer(jj).k(this.talign_t) = -1;
%                                 % note: lkhd_map is a cell, therefore we
%                                 % keep all the old cell (empty cell) but
%                                 % the length of lkhd_map will not decrease.
%                                 % For pos, z, k, they are arrays, so their
%                                 % length is constant (except the first few
%                                 % steps). This should be noticed when
%                                 % deciding the index of the element to
%                                 % be removed
%                                 this.buffer(jj).lkhd_map{this.talign_t} = [];
%                             end
%                         end
%                     end
%                 end
%                 
%                 this.talign_t = tmp_t;
%                 this.dbf_map = tmp_map;
%                 this.dbf_map = this.dbf_map/sum(sum(this.dbf_map));
%             end
%         end
%         
%         function this = updMap(this,inPara)
%             % update own map using its own measurement
%             selection = inPara.selection;
%             
%             % consensus map
%             if (selection == 1) || (selection == 3)
%                 this.cons_map = this.cons_map.*this.lkhd_map;
%             elseif (selection == 2) || (selection == 4)
%                 upd_matrix = this.upd_matrix{1};
%                 tmp_map = (this.cons_map)';
%                 tmp_map2 = upd_matrix*tmp_map(:);
%                 tmp_map = (reshape(tmp_map2,size(tmp_map)))';
%                 this.cons_map = tmp_map.*this.lkhd_map;
%             end
%             this.cons_map = this.cons_map/sum(sum(this.cons_map));
%             
%             % centralized map
%             if (selection == 1) || (selection == 3)
%                 this.cent_map = this.cent_map.*this.lkhd_map;
%             elseif (selection == 2) || (selection == 4)
%                 upd_matrix = this.upd_matrix{1};
%                 tmp_map = (this.cent_map)';
%                 tmp_map2 = upd_matrix*tmp_map(:);
%                 tmp_map = (reshape(tmp_map2,size(tmp_map)))';
%                 this.cent_map = tmp_map.*this.lkhd_map;
%             end
%             this.cent_map = this.cons_map/sum(sum(this.cent_map));
%         end
%         
%         function this = cons(this,inPara)
%             %% %%%%%%%%%%%%%%  Consensus Method %%%%%%%%%%%%%%%%%%
%             % steps:
%             % (1) observe and update the probability map for time k
%             % (2) send/receive the probability map for time k from neighbors
%             % (3) repeat step (1)            
%             
%             % update using new observation
%             %{
%             if (Selection == 1) || (Selection == 3)
%                 % update probability map
%                 for i=1:NumOfRobot
%                     tmp_cons_map = rbt_cons(i).map.*rbt_cons(i).prob;
%                     rbt_cons(i).map = tmp_cons_map/sum(sum(tmp_cons_map));
%                 end
%                 %
%             elseif (Selection == 2) || (Selection == 4)
%                 for i=1:NumOfRobot
%                     tmp_cons_map = rbt_cons(i).map;
%                     % prediction step
%                     tmp_map_cons2 = zeros(size(tmp_cons_map));
%                     for k = 1:size(pt,1)
%                         tmp_map_cons2 = tmp_map_cons2+upd_cell1{k,model_idx}*tmp_cons_map(pt(k,1),pt(k,2));
%                     end
%                     tmp_cons_map = tmp_map_cons2;
%                     
%                     % update step
%                     tmp_cons_map = tmp_cons_map.*rbt_cons(i).prob;
%                     rbt_cons(i).map = tmp_cons_map/sum(sum(tmp_cons_map));
%                 end
%             end
%             %}
%             rbt_nbhd_set = inPara.rbt_nbhd_set;
%             cons_fig = inPara.cons_fig;
%             % consensus step
%             % receive and weighted average neighboring maps            
%              
%             neigh_num=length(rbt_nbhd_set);            
%             for t = 1:neigh_num
%                 this.cons_map = this.cons_map +rbt_nbhd_set{t}.cons_map;
%             end
%             this.cons_map = this.cons_map/(neigh_num+1);
% 
%             % plot local PDFs after concensus
%             
%             if cons_fig
%                 figure
% %                 subplot(2,3,i); 
%                 contourf((this.cons_map)'); title(['Sensor ',this.idx]);
%                 % hold on;                
%                 plot(this.pos(1), this.pos(2), 's','MarkerSize',8,'LineWidth',3);       
%             end
%             
%             
% %             for i=1:NumOfRobot % Robot Iteration
% %                 rbtCon(i).map=rbt_cons(i).map;
% %             end
% %             
% %             for ConStep=1:ConsenStep % Consensus cycle
% %                 if ConsenFigure==1
% %                     fig_cnt = fig_cnt+1;
% %                     h_cons = figure(fig_cnt);
% %                     clf(h_cons);
% %                 end
% %                 for i=1:NumOfRobot % Robot Iteration
% %                     neighNum=length(this(i).neighbour)+1;
% %                     tempRbtCon(i).map = rbtCon(i).map;
% %                     for t=this(i).neighbour
% %                         tempRbtCon(i).map=tempRbtCon(i).map+rbtCon(t).map;
% %                     end
% %                     tempRbtCon(i).map=tempRbtCon(i).map/neighNum;
% %                 end
% %                 % plot local PDFs after concensus
% %                 for i=1:NumOfRobot
% %                     rbtCon(i).map=tempRbtCon(i).map;
% %                     if ConsenFigure==1
% %                         figure(fig_cnt)
% %                         subplot(2,3,i); contourf((rbtCon(i).map)'); title(['Sensor ',num2str(i)]);
% %                         hold on;
% %                         for j=1:NumOfRobot
% %                             if i==j
% %                                 plot(this(j).x, this(j).y, 's','Color',this(j).color,'MarkerSize',8,'LineWidth',3);
% %                             else
% %                                 plot(this(j).x, this(j).y, 'p','Color',this(j).color, 'MarkerSize',8,'LineWidth',1.5);
% %                             end
% %                         end
% %                     end
% %                 end
% %             end
% %             for i=1:NumOfRobot % Robot Iteration
% %                 rbt_cons(i).map=rbtCon(i).map;
% %             end
%         end
%         
%         function this = CF(this,inPara)
%             %% %%%%%%%%%%%%%% Centralized BF %%%%%%%%%%%%%%%%%%
%             % steps:
%             % (1) receive all robots' observations
%             % (2) update the probability map for time k
%             % (3) repeat step (1)
%             
%             selection = inPara.selection;
%             if (selection == 1) || (selection == 3)
%                 this.cent_map = this.cent_map.*this.lkhd_map;
%             elseif (selection == 2) || (selection == 4)
%                 % prediction step
%                 upd_matrix = this.upd_matrix{1};
%                 tmp_map = this.cent_map;
%                 tmp_map2 = upd_matrix*tmp_map(:);
%                 tmp_map = reshape(tmp_map2,size(tmp_map));
%                 
%                 % update step
%                 this.cent_map = tmp_map.*this.lkhd_map;
%             end
%             this.cent_map = this.cons_map/sum(sum(this.cent_map));            
%             
%             
% %             tmp_cent_map = rbt_cent.map;
% %             if (Selection == 1) || (Selection == 3)
% %                 % update step
% %                 for i = 1:NumOfRobot
% %                     tmp_cent_map = tmp_cent_map.*rbt_cent.prob{i};
% %                 end
% %                 rbt_cent.map = tmp_cent_map/sum(sum(tmp_cent_map));
% %                 
% %             elseif (Selection == 2) || (Selection == 4)
% %                 % prediction step
% %                 tmp_cent_map2 = zeros(size(tmp_cent_map));
% %                 for k = 1:size(pt,1)
% %                     tmp_cent_map2 = tmp_cent_map2+upd_cell1{k,model_idx}*tmp_cent_map(pt(k,1),pt(k,2));
% %                 end
% %                 tmp_cent_map = tmp_cent_map2;
% %                 
% %                 % update step
% %                 for i=1:NumOfRobot
% %                     tmp_cent_map = tmp_cent_map.*rbt_cent.prob{i};
% %                 end
% %                 rbt_cent.map = tmp_cent_map/sum(sum(tmp_cent_map));
% %             end
%         end      
%                 
%         function this = robotMove(this)
%             % needs to write up
%         end
        
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