% function [clt_res,hp_pt] = mapCluster(agent,prob_map,thresh)
function [clt_res,prob_map] = mapCluster(agent,prob_map,clt_num)
% hierarchical clustering for the probability map
% tmp_map = prob_map > thresh;
% idx = find(tmp_map);
% [x,y] = ind2sub(size(tmp_map),idx); % high probability coordinates
% hp_pt = [x,y];
% clt_res = clusterdata(hp_pt,'maxclust',2);% clustering result
% scatter3(hp_pt(:,1),hp_pt(:,2),clt_res);

% k-means clustering 
clt_idx = kmeans(prob_map(1:2,:)',clt_num);
if size(clt_idx,1) ~= 1
    % make clt_idx the row vector
    clt_idx = clt_idx';
end
% clt_idx contains the unique points and the corresponding cluster index
clt_res = [prob_map(1:2,:);clt_idx];
prob_map(4,:) = clt_idx;