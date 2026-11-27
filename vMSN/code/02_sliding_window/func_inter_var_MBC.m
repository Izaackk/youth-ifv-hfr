function variability_Inter = func_inter_var_MBC(mbc_all)
% Input: Fisher-z-transformed MSN matrices, ROI x ROI x subject.
% The self-connection is retained, matching the original analysis.
assert(size(mbc_all, 1) == size(mbc_all, 2) && size(mbc_all, 3) >= 2, ...
    'Supply square matrices for at least two subjects.');
R_ave = zeros(1, size(mbc_all, 1));
count = 0;
for i = 1:size(mbc_all, 3)
    for j = i + 1:size(mbc_all, 3)
        tmp = matcorr(mbc_all(:, :, i)', mbc_all(:, :, j)');
        tmp(isnan(tmp)) = 0;
        R_ave = R_ave + tmp;
        count = count + 1;
    end
end
variability_Inter = 1 - R_ave ./ count;
end
