function save_mgh_HFR_fs5(data, output_dir, varargin)
% Export the original HFR87-to-116-label mapping on fsaverage5 surfaces.
% An optional logical mask suppresses unselected ROI values.
template_dir = '.../HFR_output/MatchRate1/GrpTemplate_Matched_ROIs';

if nargin > 2
    data(~varargin{1}) = 0;
end
template = load(fullfile(template_dir, 'Matched87ROIs_From116ROIs.mat'), 'MatchMask');
data_87 = zeros(116, 1);
data_87(template.MatchMask) = data(:);

lh = load_mgh(fullfile(template_dir, 'Include', 'lh_Include_ROIs_FS5_lh.mgh'));
rh = load_mgh(fullfile(template_dir, 'Include', 'rh_Include_ROIs_FS5_rh.mgh'));
lh_data = zeros(numel(lh), 1);
rh_data = zeros(numel(rh), 1);
for i = 1:116
    lh_data(lh == i) = data_87(i);
    rh_data(rh == i) = data_87(i);
end

if ~isfolder(output_dir), mkdir(output_dir); end
save_mgh(lh_data, fullfile(output_dir, 'lh_fs5.mgh'), eye(4));
save_mgh(rh_data, fullfile(output_dir, 'rh_fs5.mgh'), eye(4));
end
