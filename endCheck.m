% 3/17/15
% check if the game ends

function game_end = endCheck(inPara)
prob_thresh = inPara.prob_thresh;
prob_map_pf = inPara.prob_map_pf;
r = inPara.r;
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
       