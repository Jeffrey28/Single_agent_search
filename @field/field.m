% Defines field class
% Field is assumed to be a rectange with endpoints=[xMin xMax yMin yMax]
classdef field
    
    properties
        targetNum; % specify target number
        targetPos; % specify target position, [x; y;]
        agentNum;
        step; % width of grid box in probability distribution
        endpoints; % endpoints of field, [xMin xMax yMin yMax]
        probMap;
        obs_info; % information of obstacles
    end
   
    methods
        function obj=field(endpoints,step,targetPos)
            if nargin > 0
                obj.endpoints = endpoints;
                obj.step = step;
            end
            if nargin > 2
                obj.targetPos = targetPos;
                obj.targetNum = length(targetPos(1,:));
            end
        end       
    end
end
        