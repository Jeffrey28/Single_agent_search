function [clt_res,hp_pt] = mapCluster(agent,prob_map,thresh)
% hierarchical clustering for the probability map
tmp_map = prob_map > thresh;
idx = find(tmp_map);
[x,y] = ind2sub(size(tmp_map),idx); % high probability coordinates
hp_pt = [x,y];
clt_res = clusterdata(hp_pt,'maxclust',2);% clustering result
scatter3(hp_pt(:,1),hp_pt(:,2),clt_res);