% 3/17/15
% check if the game ends

function game_end = endCheck(inPara)
prob_thresh = inPara.prob_thresh;
prob_map = inPara.prob_map;
r = inPara.r;
campus = inPara.campus;
xMin = campus.endpoints(1);
% xMax = campus.xMax;
yMin = campus.endpoints(3);
% yMax = campus.yMax;
step = campus.grid_step;

cur_cor = r.currentPos(1:2);
sig = sqrt(r.sigma_s(1,1));

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

if tot_prob > prob_thresh
    game_end = 1;
else
    game_end = 0;
end
       