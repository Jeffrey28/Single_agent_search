function [pod,pond,sim_reading] = sensorSim(agent,x,x_r) 
% simulate sensor result
% input: target position x, robot position x_r
% output: probability of detection/non-detection: pod, pond; a simulated
% senosr reading: sim_reading
sigma = agent.sigma_s;
k = agent.k_s;
% probability of detection
% x_r is robot position
pod = k/(2*pi)*(det(sigma))^(-1/2)*exp(-1/2*(x-x_r)'/sigma*(x-x_r));
pond = 1-pod;
sim_reading = (rand(1,1) < pod); % generate simulated sensor reading