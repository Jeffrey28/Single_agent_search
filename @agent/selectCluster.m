% function cur_clt = selectCluster (agent,grid_step,prob_map)
function cur_clt = selectCluster (agent,prob_map_pf)
% hp_pt = agent.hp_pt;
clt_res = prob_map_pf(4,:);
cur_clt = agent.cur_clt;

clt_num = max(agent.clt_res(3,:));
clt_idx_set = [];
% max_prob = max(prob_map(:));
max_cnt = max(prob_map_pf(3,:));

% if the current cluster still has high information, do not change cluster
if cur_clt ~= 0
    %{
    tmp_pt = hp_pt(clt_res == cur_clt,:);
    idx = sub2ind(size(prob_map),tmp_pt(:,1),tmp_pt(:,2));
    tmp_max_prob = max(prob_map(idx));
    if tmp_max_prob > max_prob/10 % may add another condition on the average probability mass
        return
    end
    %}
    tmp_cnt = prob_map_pf(3,clt_res == cur_clt);
    tmp_max_cnt = max(tmp_cnt);
    if tmp_max_cnt > max_cnt/5
        return
    end
end

% choose the nearest cluster
tmp_min_dis = []; % saves the minimum distance from the robot to each cluster
for ii = 1:clt_num
%     tmp_pt = hp_pt(clt_res == ii,:);
%     idx = sub2ind(size(prob_map),tmp_pt(:,1),tmp_pt(:,2));
%     tmp_max_prob = max(prob_map(idx));
    tmp_cnt = max(prob_map_pf(3,clt_res == ii));
    %{
    if tmp_max_prob > max_prob/10
        clt_idx_set = [clt_idx_set,ii];
%         tmp_pt = hp_pt(clt_res == ii,:);
        tmp_vec = tmp_pt*grid_step-ones(size(tmp_pt,1),1)*agent.currentPos(1:2)';
        tmp_min_dis = [tmp_min_dis,min(sum(tmp_vec.*tmp_vec,2))];
    end
    %}
    if tmp_cnt > max_cnt/5
        clt_idx_set = [clt_idx_set,ii];
        tmp_pt = prob_map_pf(1:2,clt_res == ii);
        tmp_vec = tmp_pt-agent.currentPos(1:2)*ones(1,size(tmp_pt,2));
        tmp_min_dis = [tmp_min_dis,min(sum(tmp_vec.*tmp_vec,1))];
    end
end
tmp_min_idx = find(tmp_min_dis == min(tmp_min_dis));
cur_clt = tmp_min_idx(1);