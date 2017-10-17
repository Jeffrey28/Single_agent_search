classdef Field
    properties
        fld_cor; % length of x,y coordinate  
        target; % a struct for target   
        dt; % simulation time interval
    end
    
    methods
        function this = Field(inPara)
            this.fld_cor = inPara.fld_cor;
%             this.map = ones(this.fld_size(1),this.fld_size(2))/(this.fld_size(1)*this.fld_size(2));
            this.target = inPara.target;
%             obj.tar_move = inPara.tar_move;
            this.dt = inPara.dt;
        end
        
        function this = targetMove(this)
            tar = this.target;
            f = tar.f;
            Q = tar.Q; % covariance matrix of process noise
            
            if strcmp(tar.target_model,'ped')
                % use truncated gaussian to generate theta and spd
                tmp = f(tar.state);                
                
                % here I use a very simple code, utilizing the fact that
                % noises of states are uncorrelated                
                % simple implementation
                %
                new_state = zeros(size(tmp));
                new_state(tar.pos_idx) = mvnrnd(tmp(tar.pos_idx),Q(tar.pos_idx,tar.pos_idx))';                
%                 tmp2 = normrnd(tmp(3),sqrt(Q(3,3)));
                % shift theta to be within [0,2*pi]
                pd = makedist('normal');
                pd.mu = tmp(3);
                pd.sigma = sqrt(Q(3,3));
                t1 = truncate(pd,tar.theta_bd(1),tar.theta_bd(2));
                new_state(3) = random(t1,1,1);
%                 new_state(3) = tmp2-2*pi*floor(tmp2/(2*pi));
                % use truncated gaussian to sample from range [0,inf] for
                % spd
%                 pd = makedist('normal');
                pd.mu = tmp(4);
                pd.sigma = sqrt(Q(4,4));
                t2 = truncate(pd,tar.v_bd(1),tar.v_bd(2));
                new_state(4) = random(t2,1,1);
                tar.state = new_state;
                %}
                
                % use package
                %{
%                 tmpA = [eye(4);-eye(4)];
%                 tmpB = [this.fld_cor(2);this.fld_cor(4);2*pi;3;...
%                     this.fld_cor(1);this.fld_cor(3);...
%                     0;0];
%                 tar.state = rmvnrnd(tmp,Q,1,tmpA,tmpB)';                
                l = [this.fld_cor(1);this.fld_cor(3);...
                    0;0]-tmp;
                u = [this.fld_cor(2);this.fld_cor(4);2*pi;2]-tmp;
                tmp2 = mvrandn(l,u,Q,1); % generate truncated mvn with zero mean
                tar.state = tmp2+tmp;
                %}
                tar.traj = [tar.traj,tar.state(1:2)];
            else
%                 tar.state = mvnrnd(f(tar.state),Q)';
                tar.pos = mvnrnd(f(tar.pos),Q)';
                tar.traj = [tar.traj,tar.pos];
            end
            
            
            this.target = tar;
        end
    end   
end