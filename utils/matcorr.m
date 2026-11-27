function C = matcorr(s_series, t_series)
% Correlate corresponding columns of two equally sized matrices.
% Edited by Jianxun Ren, 20171110; retained from the original workflow.
s_series = bsxfun(@minus, s_series, mean(s_series, 1));
s_series = bsxfun(@times, s_series, 1 ./ sqrt(sum(s_series.^2, 1)));
t_series = bsxfun(@minus, t_series, mean(t_series, 1));
t_series = bsxfun(@times, t_series, 1 ./ sqrt(sum(t_series.^2, 1)));
C = sum(bsxfun(@times, s_series, t_series));
end
