clear; clc;

ROIPath = '.../HFR_output/MatchRate1';
SubIDs = textread('.../SubList.txt', '%s');

script_dir = fileparts(mfilename('fullpath'));
vifp_root = fileparts(fileparts(script_dir));
output_root = fullfile(vifp_root, 'outputs', 'CmpDice');
output_all601 = fullfile(output_root, 'All_601');
if ~exist(output_all601, 'dir')
    mkdir(output_all601);
end

load(fullfile(ROIPath, 'GrpTemplate_Matched_ROIs', 'AllSelected_Patches_lh.mat'));
PatchNum_lh = AllSelected_Patches;
load(fullfile(ROIPath, 'GrpTemplate_Matched_ROIs', 'AllSelected_Patches_rh.mat'));
PatchNum_rh = AllSelected_Patches;
num_roi = sum(PatchNum_lh) + sum(PatchNum_rh);

load(fullfile(ROIPath, 'Indi_Matched_ROIs', 'All_ROI_Ind_Big.mat'));
num_sub = length(All_ROI_Ind_Big);
if num_sub ~= numel(SubIDs)
    error('Subject count mismatch: All_ROI_Ind_Big has %d, SubIDs has %d.', num_sub, numel(SubIDs));
end

% ROI indices are concatenated fsaverage4 left+right vertex indices.
num_vertices = 0;
for sub = 1:num_sub
    for roi = 1:num_roi
        idx = All_ROI_Ind_Big{sub}{roi};
        if ~isempty(idx)
            num_vertices = max(num_vertices, max(idx(:)));
        end
    end
end

dice_overlap_results = zeros(num_sub, num_sub, num_roi);
cmp_dice_results = zeros(num_sub, num_sub, num_roi);
dice_overlap_averages = zeros(num_roi, 1);
vifp_averages = zeros(num_roi, 1);
lower_tri = tril(true(num_sub), -1);

tic
for roi = 1:num_roi
    fprintf('Computing ROI %d/%d\n', roi, num_roi);

    masks = false(num_vertices, num_sub);
    for sub = 1:num_sub
        idx = All_ROI_Ind_Big{sub}{roi};
        idx = idx(~isnan(idx) & idx > 0);
        masks(idx, sub) = true;
    end

    roi_sizes = sum(masks, 1);
    intersections = double(masks') * double(masks);
    denominator = roi_sizes' + roi_sizes;

    dice_overlap = 2 * intersections ./ denominator;
    dice_overlap(denominator == 0) = NaN;
    dice_overlap(1:num_sub+1:end) = 0;

    vifp = 1 - dice_overlap;
    vifp(1:num_sub+1:end) = 0;

    dice_overlap_results(:, :, roi) = dice_overlap;
    cmp_dice_results(:, :, roi) = vifp;

    dice_values = dice_overlap(lower_tri);
    vifp_values = vifp(lower_tri);
    dice_overlap_averages(roi) = mean(dice_values(~isnan(dice_values)));
    vifp_averages(roi) = mean(vifp_values(~isnan(vifp_values)));
end
elapsed_seconds = toc;

save(fullfile(output_root, 'dice_overlap_results_new.mat'), ...
    'dice_overlap_results', 'SubIDs', 'elapsed_seconds', '-v7.3');

save(fullfile(output_root, 'cmp_dice_results_new.mat'), ...
    'cmp_dice_results', 'SubIDs', 'elapsed_seconds', '-v7.3');

% Legacy reproduction: matches the existing CmpDice/All_601/dice.mat.
averages = dice_overlap_averages;
measure = 'dice_overlap';
save(fullfile(output_all601, 'dice.mat'), ...
    'averages', 'dice_overlap_averages', 'measure', 'SubIDs', '-v7');

% Corrected vIFP/cmpdice result: vIFP = 1 - Dice overlap.
averages = vifp_averages;
measure = 'vIFP_1_minus_dice_overlap';
save(fullfile(output_all601, 'vifp.mat'), ...
    'averages', 'vifp_averages', 'dice_overlap_averages', 'measure', 'SubIDs', '-v7');

reference_file = fullfile(ROIPath, 'CmpDice', 'All_601', 'dice.mat');
if exist(reference_file, 'file')
    Ref = load(reference_file, 'averages');
    max_abs_diff = max(abs(dice_overlap_averages(:) - Ref.averages(:)));
    fprintf('Legacy dice.mat max abs diff vs reference: %.16g\n', max_abs_diff);
end

fprintf('Saved outputs to %s\n', output_root);
