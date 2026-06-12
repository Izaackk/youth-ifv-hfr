clear; clc;

% Build primary win20/step10 sliding-window vIFP.
% vIFP/cmpdice is defined as 1 - Dice overlap.

result_root = '.../HFR_output/MatchRate1';
load(fullfile(result_root, 'CmpDice', 'cmp_dice_results_new.mat'));
dice_overlap_results = cmp_dice_results;
vifp_results = 1 - dice_overlap_results;

num_sub = size(vifp_results, 1);
for roi = 1:size(vifp_results, 3)
    tmp = vifp_results(:, :, roi);
    tmp(1:num_sub+1:end) = 0;
    vifp_results(:, :, roi) = tmp;
end

SIDs = textread('.../HCPd_info/afterQC_sun.txt', '%s');
Data_info = readtable('.../HCPd_info/Subinfo_HCPD_2.xlsx');

Allsubs = Data_info.SubID;
Allage = Data_info.Age;
Allsex = Data_info.Sex;
AllmeanFD = Data_info.meanFD;

age = zeros(numel(SIDs), 1);
sex = zeros(numel(SIDs), 1);
meanFD = zeros(numel(SIDs), 1);
for nsub = 1:numel(SIDs)
    sub = SIDs{nsub};
    logi_age_sub = ismember(Allsubs, sub);
    age(nsub, 1) = Allage(logi_age_sub);
    sex(nsub, 1) = Allsex(logi_age_sub);
    meanFD(nsub, 1) = AllmeanFD(logi_age_sub);
end

[~, I] = sort(age);

window_size = 20;
step_size = 10;
num_windows = floor((numel(age) - window_size) / step_size) + 1;
num_roi = size(vifp_results, 3);

vifp_bins = zeros(num_windows, num_roi);
dice_bins = zeros(num_windows, num_roi); % kept for existing R scripts
age_bins = zeros(num_windows, 1);
sex_de_bins = zeros(num_windows, 1);
mFD_bins = zeros(num_windows, 1);
sub_win_all = cell(num_windows, 1);

lower_tri = tril(true(window_size), -1);
lower_tri_3d = repmat(lower_tri, [1, 1, num_roi]);
element_count = window_size * (window_size - 1) / 2;

for nwin = 1:num_windows
    start_idx = (nwin - 1) * step_size + 1;
    end_idx = start_idx + window_size - 1;
    I_win = I(start_idx:end_idx);
    sub_win_all{nwin} = SIDs(I_win);

    age_bins(nwin, 1) = mean(age(I_win));
    sex_de_bins(nwin, 1) = mean(sex(I_win)) - 0.5;
    mFD_bins(nwin, 1) = mean(meanFD(I_win));

    vifp_win = vifp_results(I_win, I_win, :);
    sums = sum(sum(vifp_win .* lower_tri_3d, 1), 2);
    vifp_bins(nwin, :) = squeeze(sums / element_count);
end

dice_bins = vifp_bins;
[Rmap, p] = corr(vifp_bins, age_bins);
wb_mean_vifp = mean(vifp_bins, 2);

script_dir = fileparts(mfilename('fullpath'));
vifp_root = fileparts(fileparts(script_dir));
output_dir = fullfile(vifp_root, 'repro_outputs', 'CmpDice_age');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

save(fullfile(output_dir, 'data_bin20_step10.mat'), ...
    'vifp_bins', 'dice_bins', 'age_bins', 'sex_de_bins', 'mFD_bins', ...
    'wb_mean_vifp', 'Rmap', 'p', 'window_size', 'step_size', 'sub_win_all', '-v7');

fprintf('Saved win20_step10 vIFP data to %s\n', fullfile(output_dir, 'data_bin20_step10.mat'));
