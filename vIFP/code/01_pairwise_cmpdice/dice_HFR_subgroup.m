clear; clc;

% Compute subgroup-level legacy Dice overlap and corrected vIFP.
% vIFP/cmpdice = 1 - Dice overlap.

result_root = '.../HFR_output/MatchRate1';
SubIDs = textread('.../HCPd_info/afterQC_sun.txt', '%s');

% Change this list to afterQC_sun_younger.txt if the younger subgroup is
% needed. The original script used the older subgroup.
SubIDs_subgroup = textread('.../HCPd_info/afterQC_sun_older.txt', '%s');
logi_subgroup = ismember(SubIDs, SubIDs_subgroup);

load(fullfile(result_root, 'CmpDice', 'cmp_dice_results_new.mat'), 'cmp_dice_results');
dice_overlap_results = cmp_dice_results(logi_subgroup, logi_subgroup, :);

num_sub = numel(SubIDs_subgroup);
num_roi = size(dice_overlap_results, 3);
lower_tri = tril(true(num_sub), -1);

dice_overlap_averages = zeros(num_roi, 1);
vifp_averages = zeros(num_roi, 1);
for roi = 1:num_roi
    dice_overlap = dice_overlap_results(:, :, roi);
    dice_values = dice_overlap(lower_tri);
    dice_values = dice_values(~isnan(dice_values));

    dice_overlap_averages(roi) = mean(dice_values);
    vifp_averages(roi) = mean(1 - dice_values);
end

averages = vifp_averages;
