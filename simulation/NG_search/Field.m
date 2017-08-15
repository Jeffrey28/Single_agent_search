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
            tar.pos = f(tar.pos)+mvnrnd([0,0],Q)';
            
            % obsolete
%             A = tar.A; % A matrix
%             B = tar.B;            
%             tar.pos = A*tar.pos+B+mvnrnd([0,0],Q)';

            tar.traj = [tar.traj,tar.pos];
            this.target = tar;
        end
    end   
end