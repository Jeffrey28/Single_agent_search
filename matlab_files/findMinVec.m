function tmp_idx = findMinVec(mat,vec)
% find the indices of closest vectors in matrix mat to vec
n = size(mat,2);
tmp_dif_vec = mat-vec*ones(1,n);
tmp_sum = sqrt(sum(tmp_dif_vec.*tmp_dif_vec,1));
tmp_idx = tmp_sum == min(tmp_sum);
    