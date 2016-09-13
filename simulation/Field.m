classdef Field
    properties
        map; % field map
        fld_size; % length of x,y coordinate  
        target; % a struct for target   
%         tar_move; % indicator whether the target moves
        dt; % simulation time interval
    end
    
    methods
        function this = Field(inPara)
            this.fld_size = inPara.fld_size;
%             this.map = ones(this.fld_size(1),this.fld_size(2))/(this.fld_size(1)*this.fld_size(2));
            this.target = inPara.target;
%             obj.tar_move = inPara.tar_move;
            this.dt = inPara.dt;
        end
        
        function this = targetMove(this)
            tar = this.target;
            A = tar.A; % A matrix
            Q = tar.Q; % covariance matrix of process noise
            tar.pos = A*tar.pos+mvnrnd([0,0],Q);
            this.target = tar;
        end
            
%         function fld = setMotionModel(fld)
%             % later will change to linear state space model
%             if fld.tar_move == 0
%                 fld.target.dx= 0;
%                 fld.target.dy= 0;
%             elseif fld.tar_move == 1
%                 model_idx = fld.target.model_idx;
%                 fld.target.dx= fld.target.dx_set(model_idx);
%                 fld.target.dy= fld.target.dx_set(model_idx);
%             end         
%         end
        

%         function this = targetMove(this)
%             
%             tmp_u_set = this.target.u_set;
%             tmp_idx = this.target.model_idx;
%             if this.tar_move == 1
%                 % Target Moves
%                 % if current model makes target out of field, choose the next
%                 % motion model in fld.target.dx_set and fld.target.dy_set
%                 tmp_pos = this.target.pos + tmp_u_set(tmp_idx)*this.dt;                
%                 while (tmp_pos(1) <= 0 || tmp_pos(2) <= 0 || tmp_pos(1) >= this.fld_size(1) || tmp_pos(2) >= this.fld_size(2))
%                     tmp_idx = rem(tmp_idx+1,length(this.target.u_set));
%                     if tmp_idx == 0
%                         tmp_idx = length(this.target.u_set);
%                     end
%                     tmp_pos = this.target.pos + tmp_u_set(tmp_idx)*this.dt;
%                 end
%                 this.target.pos = tmp_pos;
%                 this.target.traj = [this.target.traj;tmp_pos];
%                 this.target.model_idx = tmp_idx;
%                 
% %                 tmp_tx = fld.target.pos(1) + fld.target.speed * fld.target.dx;
% %                 tmp_ty = fld.target.pos(2) + fld.target.speed * fld.target.dy;
% %                 tmp_model_idx = fld.target.model_idx;
% %                 while (tmp_tx <= 1 || tmp_ty <= 1 || tmp_tx >= fld.x || tmp_ty >= fld.y)
% %                     tmp_model_idx = rem(tmp_model_idx+1,size(fld.target.dx_set,2));
% %                     if tmp_model_idx == 0
% %                         tmp_model_idx = size(fld.target.dx_set,2);
% %                     end
% %                     tmp_tx = fld.tx + fld.target.speed * fld.target.dx_set(tmp_model_idx);
% %                     tmp_ty = fld.ty + fld.target.speed * fld.target.dy_set(tmp_model_idx);
% %                 end
% %                 fld.target.pos = [tmp_tx;tmp_ty];
% %                 fld.target.model_idx = tmp_model_idx;
%             end
%         end
    end   
end