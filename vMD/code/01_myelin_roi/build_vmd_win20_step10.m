clear; clc;

% Compute ROI-wise myelin-distribution variability on fsaverage4 surfaces.
myelin_path = '.../HCP-D_surf_myelin';
parc_path = '.../HFR_output/MatchRate1/Indi_HFR_parcellataion';
mask_file = '.../templates/ind_fs4_surf_HCP.mat';
sid_file = '.../SubList.txt';
info_file = '.../SubjectInfo.xlsx';
out_root = '.../HFR_output/MatchRate1/Myelin';
window_size = 20;
step_size = 10;

validateattributes(window_size, {'numeric'}, {'scalar', 'integer', 'finite', '>=', 2});
validateattributes(step_size, {'numeric'}, {'scalar', 'integer', 'finite', 'positive'});
addpath(fileparts(mfilename('fullpath')));

% Match covariates to the supplied subject-list order, including age ties.
SIDs = strtrim(string(splitlines(fileread(sid_file))));
SIDs(SIDs == "") = [];
assert(~isempty(SIDs) && numel(unique(SIDs)) == numel(SIDs), ...
    'The subject list must be nonempty and contain unique IDs.');
info = readtable(info_file, 'TextType', 'string');
required = {'SubID', 'Age', 'Sex', 'meanFD'};
assert(all(ismember(required, info.Properties.VariableNames)), ...
    'The covariate table must contain SubID, Age, Sex, and meanFD.');
table_ids = string(info.SubID);
[found, index] = ismember(SIDs, table_ids);
assert(all(found), 'Some subjects are missing from the covariate table.');
assert(sum(ismember(table_ids, SIDs)) == numel(SIDs), ...
    'Selected subjects must have exactly one row in the covariate table.');
age = info.Age(index);
sex = info.Sex(index);
mFD = info.meanFD(index);
assert(isnumeric(age) && isnumeric(sex) && isnumeric(mFD), ...
    'Age, Sex, and meanFD must be numeric.');
assert(all(isfinite([age; sex; mFD])), 'Covariates must not contain missing values.');
assert(numel(SIDs) >= window_size, 'The subject count is smaller than the window.');
mask = load(mask_file, 'ind_mask');
ind_mask = double(mask.ind_mask(:));
assert(~isempty(ind_mask) && all(isfinite(ind_mask)) && ...
    all(ind_mask >= 1 & ind_mask <= 5124 & ind_mask == fix(ind_mask)) && ...
    numel(unique(ind_mask)) == numel(ind_mask), ...
    'ind_mask must contain unique one-based indices into [left; right] fsaverage4 vertices.');
myelin_all = zeros(numel(ind_mask), numel(SIDs));
parc_all = zeros(size(myelin_all));
for i = 1:numel(SIDs)
    sub = char(SIDs(i));
    left = gifti(fullfile(myelin_path, sub, [sub, '.L.MyelinMap_BC_MSMAll.32k_fs_LR_faverage4.func.gii']));
    right = gifti(fullfile(myelin_path, sub, [sub, '.R.MyelinMap_BC_MSMAll.32k_fs_LR_faverage4.func.gii']));
    left_parc = load_mgh(fullfile(parc_path, sub, 'lh_HFR_87.mgh'));
    right_parc = load_mgh(fullfile(parc_path, sub, 'rh_HFR_87.mgh'));
    assert(numel(left.cdata) == 2562 && numel(right.cdata) == 2562 && ...
        numel(left_parc) == 2562 && numel(right_parc) == 2562, ...
        'Each hemisphere must contain 2562 fsaverage4 vertices.');
    values = [double(left.cdata(:)); double(right.cdata(:))];
    parcels = [double(left_parc(:)); double(right_parc(:))];
    myelin_all(:, i) = values(ind_mask);
    parc_all(:, i) = parcels(ind_mask);
end
assert(~any(isinf(myelin_all(:))), 'Myelin input must not contain infinite values.');
assert(all(ismember(parc_all(:), 0:87)), 'Parcellation labels must be integers from 0 to 87.');
[~, order] = sort(age);
nwin = floor((numel(SIDs) - window_size) / step_size) + 1;
H = zeros(nwin, 87);
age_bins = zeros(nwin, 1);
sex_de_bins = zeros(nwin, 1);
mFD_bins = zeros(nwin, 1);
sub_win_all = cell(nwin, window_size);
for w = 1:nwin
    index = order((w - 1) * step_size + (1:window_size));
    H(w, :) = func_vmd_roi(myelin_all(:, index), parc_all(:, index), 87);
    age_bins(w) = mean(age(index));
    sex_de_bins(w) = mean(sex(index));
    mFD_bins(w) = mean(mFD(index));
    sub_win_all(w, :) = reshape(cellstr(SIDs(index)), 1, []);
end
vMD_bins = H;
SIDs_used = cellstr(SIDs);
roi_out = fullfile(out_root, 'ROI');
wb_out = fullfile(out_root, 'WB');
if ~exist(roi_out, 'dir'), mkdir(roi_out); end
if ~exist(wb_out, 'dir'), mkdir(wb_out); end
filename = sprintf('data_bin%d_step%d.mat', window_size, step_size);
save(fullfile(roi_out, filename), 'H', 'vMD_bins', 'age_bins', 'sex_de_bins', ...
    'mFD_bins', 'window_size', 'step_size', 'sub_win_all', 'SIDs_used', '-v7');
% Whole-brain output is the unweighted mean of the 87 ROI H statistics.
H = mean(H, 2);
save(fullfile(wb_out, filename), 'H', 'age_bins', 'sex_de_bins', ...
    'mFD_bins', 'window_size', 'step_size', 'sub_win_all', 'SIDs_used', '-v7');
fprintf('Saved vMD: %d windows x 87 ROIs, plus the ROI mean.\n', nwin);
