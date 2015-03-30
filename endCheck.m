% 3/17/15
% check if the game ends

function game_end = endCheck(inPara)
prob_thresh = inPara.prob_thresh;
prob_map_pf = inPara.prob_map_pf;
r = inPara.r;
sensor_reading = inPara.sensor_reading;
particles = inPara.particles;
% particles = inPara.particles;
% campus = inPara.campus;
% xMin = campus.endpoints(1);
% % xMax = campus.xMax;
% yMin = campus.endpoints(3);
% % yMax = campus.yMax;
% step = campus.grid_step;

cur_cor = r.currentPos(1:2);
sig = sqrt(r.sigma_s(1,1));

% this part is the version when the probability map is continuous, not
% prarticles
%{
tot_prob = 0;

for x = cur_cor(1)-sig:step:cur_cor(1)+sig
    for y = cur_cor(2)-sig:step:cur_cor(2)+sig   
        if norm(cur_cor-[x;y]) <= 0.8*sig
            % change the actual coordinates into the grid index of the
            % prob_map
            tmp_cor = [floor((x-xMin)/step)+1,floor((y-yMin)/step)+1];
            tot_prob = tot_prob+prob_map(tmp_cor(1),tmp_cor(2));
        end
    end
end
%}

% this part counts then number of particles inside the range of view
%{
tot_cnt = 0;
for ii = 1:size(prob_map_pf,2)
    if norm(cur_cor-prob_map_pf(1:2,ii)) <= 0.8*sig
       tot_cnt = prob_map_pf(3,ii);
    end
end
n_data = sum(prob_map_pf(3,:));
tot_prob = tot_cnt/n_data;
if tot_prob > prob_thresh
    game_end = 1;
else
    game_end = 0;
end     
%}

% use two conditions, when either is satisfied, the game ends

% condition 1: when the sum of densities inside the range of view is
% greater than a threshold

tmp_pt = prob_map_pf(1:2,:);
tmp_vec = tmp_pt-cur_cor*ones(1,size(prob_map_pf,2));
% find the indices of particles whose distance to the robot is less than
% 0.8*sig
%
tmp_idx = [sum(tmp_vec.*tmp_vec,1) <= sig^2*ones(1,size(tmp_pt,2))];
all_pdf_sum = sum(prob_map_pf(5,:));
fov_pdf_sum = sum(prob_map_pf(5,tmp_idx));
if fov_pdf_sum > prob_thresh*all_pdf_sum
    game_end = 1;
    return
else
    game_end = 0;
    return
end  
%}

% condition 2: 5 consecutive positive observations
%{
k = length(sensor_reading);% sensor_reading length is equal to current time
if k >= 2
    if (sum(sensor_reading(k-1:k)) == 5)
        game_end = 1;
    end
else
    game_end = 0;
end
%}
% condition 3: the covariance is below a threshold

uni_par = (unique(particles','rows'))'; % 2-by-x matrix
mean_p = mean(uni_par,2);
n_uniq = size(uni_par,2); % # of unique points
%{
mean_p = mean(particles,2);
n = size(particles,2);
dif = particles - mean_p*ones(1,n);
cov_p = dif*dif'/n;
if norm(cov_p,2) < 1
    game_end = 1;
else
    game_end = 0;
end
%}
dif = uni_par - mean_p*ones(1,n_uniq);
cov_p = dif*dif'/n_uniq;
if norm(cov_p,2) < 5
    game_end = 1;
    return
else
    game_end = 0;
    return
end