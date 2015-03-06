% Defines agent class
% Agent type describes the sensor model

classdef agent
    properties
        type; % agent type: 'human', 'robot'
        maxV; % max travel velocity: h=1.5, r=2
%         maxA; % max acceleration
        a_lb; % min deceleration
        a_ub; % max acceleration
        w_lb; % min angular velocity
        w_ub; % max angular velocity
%         dV; % incremental velocity change: h=0.5, r=2
%         maxDPsi; % max turning angle: h=pi, r=pi/2
%         dPsi; % Incremental heading change: h=pi/8, r=pi/4
%         visionR; % h=1.5, r=2
%         visionTheta; % h=pi,r=2*pi
        %probDetTrueReadingInFOV; % h=0.85, r=0.7
        %probDetFalseReadingInFOV; % h=0.1, r=0.2
        %probDetNotInFOV; % h=0.05, r=0.1
%         ProbDGivenT;
%         ProbDOutsideFOV;
%         planHorizon; % h=2,r=2
        currentPos; % current position for the agent, [x y psi]
        currentV; % current velocity
        traj; % current and planned trajectory for agent at each time step
        sigma_s; % covariance for sensor model
        k_s; % the scaling constant for sensor model
        intgr_step = 0.5; % step size for integral over x when calculating POD
%         FOVcoords; % grid coordinates contained in agent field of view
%         human
        % No entry zone: 
        % - doesnt need to be initialized
        % - noEntryZone(1,:,:) = [ array x coords; array y coords] of 
        %   polygon that contains no enty zone.
%         noEntryZone; % contains x and y corners of polygon 
        
        % Weights for cost function
%         wQz;
%         wDz;
    end
        
    methods
        
        % Contruction function
        function obj=agent(type)
            obj;
            if nargin > 0
                obj.type = lower(type);
                if nargin > 1
                    % add more terms that can be used to define agents here
                end
            end
        end
        
        % Move with bestAction
        function obj = takeNextAction(obj, bestAction)
            obj.currentPos = obj.currentPos + bestAction;
        end
        
        % More methods are defined in the agent folder, including:
        % legalActions = returnLegalActions(agent,field)
        % [bestAction,plannedSteps,plannedProb,FOV,plannedQz] = getNextAction(agent,field,nextSteps,priorProb,k,kh)
        % generateCostFunction
        % getNextActionWithDanger
        % inFieldOfView
        % sensorModel
        % updateProbability
        
    end
end