clear; clc;

% Compute age-sorted, ROI-wise MSN variability; save GAM-ready matrices.
data_path = '.../HFR_output/MatchRate1/MBC/data_HFR';
sid_file = '.../SubList.txt';
info_file = '.../SubjectInfo.xlsx';
out_root = '.../HFR_output/MatchRate1/MBC/diff_MBC_minmaxnorm_novetices/Inter';
window_size = 20;
step_size = 10;

validateattributes(window_size, {'numeric'}, {'scalar', 'integer', 'finite', '>=', 2});
validateattributes(step_size, {'numeric'}, {'scalar', 'integer', 'finite', 'positive'});
module_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
repo_root = fileparts(module_root);
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(repo_root, 'utils'));

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
MSN_all = zeros(87, 87, numel(SIDs));
for i = 1:numel(SIDs)
    loaded = load(fullfile(data_path, char(SIDs(i)), 'MBC_minmaxnorm_novertices.mat'), 'MBC');
    assert(isequal(size(loaded.MBC), [87, 87]), 'Each MSN must be 87 x 87.');
    MSN_all(:, :, i) = Fisherz(loaded.MBC);
end
[~, order] = sort(age);
nwin = floor((numel(SIDs) - window_size) / step_size) + 1;
var_all = zeros(nwin, 87);
age_bins = zeros(nwin, 1);
sex_bins = zeros(nwin, 1);
mFD_bins = zeros(nwin, 1);
sub_win_all = cell(nwin, window_size);
for w = 1:nwin
    index = order((w - 1) * step_size + (1:window_size));
    var_all(w, :) = func_inter_var_MBC(MSN_all(:, :, index));
    age_bins(w) = mean(age(index));
    sex_bins(w) = mean(sex(index));
    mFD_bins(w) = mean(mFD(index));
    sub_win_all(w, :) = reshape(cellstr(SIDs(index)), 1, []);
end
mvar_all = mean(var_all, 2);
% Aliases allow the existing ROI GAM template to consume this output directly.
H = var_all;
sex_de_bins = sex_bins;
SIDs_used = cellstr(SIDs);
if ~exist(out_root, 'dir'), mkdir(out_root); end
filename = fullfile(out_root, sprintf('data_bin%d_step%d.mat', window_size, step_size));
save(filename, 'var_all', 'mvar_all', 'H', 'age_bins', 'sex_bins', ...
    'sex_de_bins', 'mFD_bins', 'window_size', 'step_size', 'sub_win_all', 'SIDs_used', '-v7');
fprintf('Saved vMSN: %d windows x 87 ROIs.\n', nwin);
