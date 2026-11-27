clear; clc;

% Compute age-sorted, ROI-wise structural-connectome variability.
data_path = '.../HCPd_SC_HFR88indi_20260102';
sid_file = '.../SubList.txt';
info_file = '.../SubjectInfo.xlsx';
out_root = '.../HFR_output/MatchRate1/SC_DWI_20260102/Inter_new';
window_size = 20;
step_size = 10;
exclude_sids = {}; % Optional IDs to omit from the supplied subject list.

script_dir = fileparts(mfilename('fullpath'));
vsc_root = fileparts(fileparts(script_dir));
addpath(fullfile(vsc_root, 'code', 'utils'));

validateattributes(window_size, {'numeric'}, {'scalar', 'integer', 'finite', '>=', 2});
validateattributes(step_size, {'numeric'}, {'scalar', 'integer', 'finite', 'positive'});
SIDs = strtrim(string(splitlines(fileread(sid_file))));
SIDs(SIDs == "" | ismember(SIDs, string(exclude_sids))) = [];
assert(numel(unique(SIDs)) == numel(SIDs), 'Subject IDs must be unique.');
SIDs_used = cellstr(SIDs);
data_info = readtable(info_file, 'TextType', 'string');

nsub = numel(SIDs);
assert(nsub >= window_size, 'The subject count is smaller than the window.');
table_ids = string(data_info.SubID);
[found, index] = ismember(SIDs, table_ids);
assert(all(found) && sum(ismember(table_ids, SIDs)) == nsub, ...
    'Each selected subject must have exactly one row in the covariate table.');
age = data_info.Age(index);
sex = data_info.Sex(index);
mFD = data_info.meanFD(index);
assert(isnumeric(age) && isnumeric(sex) && isnumeric(mFD) && ...
    all(isfinite([age; sex; mFD])), 'Age, Sex, and meanFD must be finite numeric values.');

logi_87 = true(88, 1);
logi_87(1) = false;

SC_all = zeros(nsub, 87, 87);
for i = 1:nsub
    sub = char(SIDs(i));
    sub_dir = fullfile(data_path, sub);
    mat_file = fullfile(sub_dir, [sub, '_HFR_88_indi_SC.mat']);
    csv_file = fullfile(sub_dir, [sub, '_HFR_88_indi_SC.csv']);

    if isfile(mat_file)
        loaded = load(mat_file, 'SC');
        SC = loaded.SC;
    else
        SC = readmatrix(csv_file);
    end
    assert(isequal(size(SC), [88, 88]), 'Each input connectome must be 88 x 88.');

    SC(1:88+1:88^2) = 0;
    SC_all(i, :, :) = SC(logi_87, logi_87);
end

[~, I] = sort(age);
num_windows = floor((nsub - window_size) / step_size) + 1;

var_bins = zeros(num_windows, 87);
age_bins = zeros(num_windows, 1);
sex_bins = zeros(num_windows, 1);
mFD_bins = zeros(num_windows, 1);
wb_mean = zeros(num_windows, 1);
sub_win_all = cell(num_windows, window_size);

for nwin = 1:num_windows
    fprintf('Window %d / %d\n', nwin, num_windows);

    start_idx = (nwin - 1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    I_win = I(start_idx:end_idx);
    sub_win_all(nwin, :) = reshape(cellstr(SIDs(I_win)), 1, []);

    age_bins(nwin, 1) = mean(age(I_win));
    sex_bins(nwin, 1) = mean(sex(I_win));
    mFD_bins(nwin, 1) = mean(mFD(I_win));

    count = 0;
    ave_rmat = zeros(1, 87);

    for i = 1:numel(I_win)
        SC_i = squeeze(SC_all(I_win(i), :, :));
        for j = i + 1:numel(I_win)
            SC_j = squeeze(SC_all(I_win(j), :, :));
            tmp = matcorr(SC_i', SC_j');
            tmp(isnan(tmp)) = 0;
            ave_rmat = ave_rmat + tmp;
            count = count + 1;
        end
    end

    inter_similarity = ave_rmat ./ count;
    var_bins(nwin, :) = 1 - inter_similarity;
    wb_mean(nwin, 1) = mean(var_bins(nwin, :));
end

if ~isfolder(out_root), mkdir(out_root); end
out_file = fullfile(out_root, sprintf('data_bin%d_step%d.mat', window_size, step_size));
save(out_file, ...
    'var_bins', 'age_bins', 'sex_bins', 'mFD_bins', ...
    'wb_mean', 'window_size', 'step_size', ...
    'sub_win_all', 'SIDs_used', 'sid_file', 'exclude_sids', '-v7');

fprintf('Saved vSC: %d windows x 87 ROIs to %s\n', num_windows, out_file);
