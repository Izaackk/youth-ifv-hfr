clc;

script_dir = fileparts(mfilename('fullpath'));
vsc_root = fileparts(fileparts(script_dir));
repo_root = fileparts(vsc_root);
addpath(fullfile(vsc_root, 'code', 'utils'));

if ~exist('data_path', 'var')
    data_path = '.../HCPd_SC_HFR88indi_20260102';
end
if ~exist('info_path', 'var')
    info_path = '.../HCPd_info';
end
if ~exist('sid_file', 'var')
    sid_file = fullfile(info_path, 'afterQC_sun_afterdMRI.txt');
end
if ~exist('out_root', 'var')
    out_root = fullfile(repo_root, 'repro_outputs', 'vSC', 'SC_DWI_20260102', 'Inter_new');
end
if ~exist('window_size', 'var')
    window_size = 20;
end
if ~exist('step_size', 'var')
    step_size = 10;
end
if ~exist('out_file', 'var')
    out_file = fullfile(out_root, sprintf('data_bin%d_step%d.mat', window_size, step_size));
end

if ~exist(out_root, 'dir')
    mkdir(out_root);
end

if ~exist('exclude_sids', 'var')
    exclude_sids = {};
elseif ischar(exclude_sids)
    exclude_sids = {exclude_sids};
elseif isstring(exclude_sids)
    exclude_sids = cellstr(exclude_sids);
end

SIDs = textread(sid_file, '%s');
if ~isempty(exclude_sids)
    keep_mask = ~ismember(SIDs, exclude_sids);
    fprintf('Excluded %d subject(s) from %s\n', sum(~keep_mask), sid_file);
    SIDs = SIDs(keep_mask);
end
SIDs_used = SIDs;
data_info = readtable(fullfile(info_path, 'Subinfo_HCPD_2.xlsx'));

allsubs = data_info.SubID;
allage = data_info.Age;
allsex = data_info.Sex;
allmFD = data_info.meanFD;

nsub = numel(SIDs);
age = zeros(nsub, 1);
sex = zeros(nsub, 1);
mFD = zeros(nsub, 1);

for i = 1:nsub
    sub = SIDs{i};
    idx = ismember(allsubs, sub);
    age(i, 1) = allage(idx);
    sex(i, 1) = allsex(idx);
    mFD(i, 1) = allmFD(idx);
end

logi_87 = true(88, 1);
logi_87(1) = false;

SC_all = zeros(nsub, 87, 87);
for i = 1:nsub
    sub = SIDs{i};
    sub_dir = fullfile(data_path, sub);
    mat_file = fullfile(sub_dir, [sub, '_HFR_88_indi_SC.mat']);
    csv_file = fullfile(sub_dir, [sub, '_HFR_88_indi_SC.csv']);

    if exist(mat_file, 'file')
        load(mat_file, 'SC');
    elseif exist(csv_file, 'file')
        SC = readmatrix(csv_file);
        save(mat_file, 'SC');
    else
        error('SC file not found for %s', sub);
    end

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
    sub_win_all(nwin, :) = reshape(SIDs(I_win), 1, []);

    age_bins(nwin, 1) = mean(age(I_win));
    sex_bins(nwin, 1) = mean(sex(I_win));
    mFD_bins(nwin, 1) = mean(mFD(I_win));

    count = 0;
    ave_rmat = zeros(1, 87);

    for i = 1:numel(I_win)
        SC_i = squeeze(SC_all(I_win(i), :, :));
        for j = i + 1:numel(I_win)
            SC_j = squeeze(SC_all(I_win(j), :, :));
            tmp = my_matcorr(SC_i', SC_j');
            tmp(isnan(tmp)) = 0;
            ave_rmat = ave_rmat + tmp;
            count = count + 1;
        end
    end

    inter_similarity = ave_rmat ./ count;
    var_bins(nwin, :) = 1 - inter_similarity;
    wb_mean(nwin, 1) = mean(var_bins(nwin, :));
end

save(out_file, ...
    'var_bins', 'age_bins', 'sex_bins', 'mFD_bins', ...
    'wb_mean', 'window_size', 'step_size', ...
    'sub_win_all', 'SIDs_used', 'sid_file', 'exclude_sids', '-v7');

fprintf('Saved Step 5 vSC data to %s\n', out_file);
