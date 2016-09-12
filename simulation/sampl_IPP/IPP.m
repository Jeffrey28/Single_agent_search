% main function for running the igRRT code

clear
addpath ('common');

% generate map
range = [0;10;0;10]; % map corner coordinates: [xm;xM;ym;yM]
grid_res = 0.5; % grid resolution
tar_pos = [8;8];

% prob map
prior_map = zeros((range(2)-range(1))/grid_res,(range(2)-range(1))/grid_res);
norm_mean = [(range(1)+range(2))/2,(range(3)+range(4))/2];
norm_cov = ((range(1)+range(2))/2)^2*eye(2);
for ii = 1:size(prior_map,1)
    for jj = 1:size(prior_map,2)
        prior_map(ii,jj) = mvnpdf([ii,jj],norm_mean,norm_cov);
    end
end

% obstacle map
obs_map = {[4,4,5,5;5,9,9,5],[4,4,5,5;4,5,5,4],[5,5,9,9;4,5,5,4]};

inPara = struct('probmap', {prior_map},'verbose', false,'obsmap',{obs_map},'tarpos',tar_pos);

map = Map(inPara);

% robot model
veh = Vehicle([], 'stlim', 1.2);

% RRT for IPP
rrt = igRRT(map,veh,'range',range);
rrt.plan;
rrt.plot();
% rrt.path([1;2;0],[tar_pos;0]);

%% missing parts:
% how to use information gain to guide the search (ideally, the path will tend towards high prob region)
% some bugs exist: 
% the extended graph concentrates near origin (and there are points in negative coordinate!!)
% the best next feasible point xnew seems not correct. check if there's any
% bug here.

% how to estimate info gain
% test different probability map, have the demo

% later, think about how to use finite horizon prediction and real-time
% rewiring for better search path planning