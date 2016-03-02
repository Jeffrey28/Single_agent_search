%% convert SLAM data file .pgm to .mat file
clear all

cur_folder_path = pwd;
folder_path = strcat(cur_folder_path,'/SLAM_map');
file_name = 'laserdata_021816_map.pgm';

occu_map = imread(fullfile(folder_path,file_name));

% find the unoccupied and occupied areas
unocc_val = 254; 
occ_val = 0;
prob_map = zeros(size(occu_map));
% marker for unoccupied positions
unocc_marker = (unocc_val - occu_map) <0.01;
% marker for occupied positions
occ_marker = (occu_map-occ_val) <0.01;
% marker for positions out-of-range
unused_marker = ones(size(occu_map));
unused_marker (unocc_marker | occ_marker) = 0;
unused_marker = logical(unused_marker);

prob_map(unocc_marker) = 1;
prob_map(occ_marker) = 0;
prob_map = prob_map/sum(sum(prob_map)); % normalize probability map
prob_map(unused_marker) = -100; % assign very low information gain for out-of-range regions

map = struct();
map.occu_map = int16(occu_map);
map.prob_map = prob_map;
% find coordinates of unoccupied, occupied and out-of-range positions in
% grid map
[tmpx1,tmpy1] = find(unocc_marker);
map.unocc_cor = [tmpx1,tmpy1]';
[tmpx2,tmpy2] = find(occ_marker);
map.occ_cor = [tmpx2,tmpy2]';
[tmpx3,tmpy3] = find(unused_marker);
map.unused_cor = [tmpx3,tmpy3]';
save(fullfile(folder_path,'shared_room_map.mat'),'map')