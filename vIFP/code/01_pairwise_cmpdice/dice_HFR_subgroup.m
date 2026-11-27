clear; clc;

% Compute subgroup-level legacy Dice overlap and corrected vIFP.
% vIFP/cmpdice = 1 - Dice overlap.

pairwise_file = '.../vIFP/outputs/CmpDice/cmp_dice_results_new.mat';
subgroup_file = '.../SubList_subgroup.txt';

Pairwise = load(pairwise_file);
assert(isfield(Pairwise, 'measure') && strcmp(Pairwise.measure, 'vIFP_1_minus_dice_overlap'), ...
    'Use the marked vIFP output from dice_HFR_group, not a legacy Dice-overlap file.');
SubIDs = cellstr(string(Pairwise.SubIDs(:)));
SubIDs_subgroup = cellstr(strtrim(string(splitlines(fileread(subgroup_file)))));
SubIDs_subgroup(cellfun(@isempty, SubIDs_subgroup)) = [];
assert(numel(unique(SubIDs_subgroup)) == numel(SubIDs_subgroup) && ...
    all(ismember(SubIDs_subgroup, SubIDs)), 'Subgroup IDs must be unique and present in the pairwise file.');
logi_subgroup = ismember(SubIDs, SubIDs_subgroup);

vifp_results = Pairwise.cmp_dice_results(logi_subgroup, logi_subgroup, :);

num_sub = sum(logi_subgroup);
assert(num_sub >= 2, 'A subgroup must contain at least two subjects.');
num_roi = size(vifp_results, 3);
lower_tri = tril(true(num_sub), -1);

dice_overlap_averages = zeros(num_roi, 1);
vifp_averages = zeros(num_roi, 1);
for roi = 1:num_roi
    vifp = vifp_results(:, :, roi);
    vifp_values = vifp(lower_tri);
    vifp_values = vifp_values(~isnan(vifp_values));

    vifp_averages(roi) = mean(vifp_values);
    dice_overlap_averages(roi) = mean(1 - vifp_values);
end

averages = vifp_averages;
