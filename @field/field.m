% Defines field class
% Field is assumed to be a rectange with endpoints=[xMin xMax yMin yMax]
classdef field
    
    properties
        targetNum; % specify target number
        targetPos; % specify target position, [x; y;]
        agentNum;
        grid_step; % width of grid box in probability distribution
        endpoints; % endpoints of field, [xMin xMax yMin yMax]
        obs_info; % information of obstacles
        w; % mixture weight for gaussian mixture model
        mu; % vector of mean
        sigma; % matrix of covariance 
    end
   
    methods
        function obj=field(endpoints,step,targetPos,w,mu,sigma) 
            obj.endpoints = endpoints;
            obj.grid_step = step;      
            obj.targetPos = targetPos;
            obj.targetNum = length(targetPos(1,:));
            obj.w = w;
            obj.mu = mu;
            obj.sigma = sigma;
        end       
    end
end
        