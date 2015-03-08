function cur_clt = selectCluster (agent,grid_step,prob_map)
hp_pt = agent.hp_pt;
clt_res = agent.clt_res;
cur_clt = agent.cur_clt;

clt_num = max(agent.clt_res);
clt_idx_set = [];
max_prob = max(prob_map(:));

% if the current cluster still has high information, do not change cluster
if cur_clt ~= 0
    tmp_pt = hp_pt(clt_res == cur_clt,:);
    idx = sub2ind(size(prob_map),tmp_pt(:,1),tmp_pt(:,2));
    tmp_max_prob = max(prob_map(idx));
    if tmp_max_prob > max_prob/10 % may add another condition on the average probability mass
        return
    end
end

% choose the nearest cluster
tmp_min_dis = [];
for ii = 1:clt_num
    tmp_pt = hp_pt(clt_res == ii,:);
    idx = sub2ind(size(prob_map),tmp_pt(:,1),tmp_pt(:,2));
    tmp_max_prob = max(prob_map(idx));
    if tmp_max_prob > max_prob/10
        clt_idx_set = [clt_idx_set,ii];
%         tmp_pt = hp_pt(clt_res == ii,:);
        tmp_vec = tmp_pt*grid_step-ones(size(tmp_pt,1),1)*agent.currentPos(1:2)';
        tmp_min_dis = [tmp_min_dis,min(sum(tmp_vec.*tmp_vec,2))];
    end
end
tmp_min_idx = find(tmp_min_dis == min(tmp_min_dis));
cur_clt = tmp_min_idx;