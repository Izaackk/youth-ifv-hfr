clear; clc;

% Extract eight morphometric features in the manuscript's HFR87 order.
parc_path = '.../HFR_output/MatchRate1/Indi_HFR_parcellataion';
sid_file = '.../SubList.txt';
out_root = '.../HFR_output/MatchRate1/MBC/data_HFR';

addpath(fileparts(mfilename('fullpath')));
SIDs = strtrim(string(splitlines(fileread(sid_file))));
SIDs(SIDs == "") = [];
fields = {'vertices', 'surface_area', 'gray_matter_vol', 'thickness_mean', ...
    'thickness_std', 'mean_curv', 'gauss_curv', 'folding_index', 'curv_index'};
% Native stats use Network_2 through Network_88 for analysis ROIs 1:87.
roi_names = cellstr(compose('Network_%d', (2:88)'));
for n = 1:numel(SIDs)
    sub = char(SIDs(n));
    [left, left_names] = read_morph_stats(fullfile(parc_path, sub, 'lh.aparc.HFR.stats'));
    [right, right_names] = read_morph_stats(fullfile(parc_path, sub, 'rh.aparc.HFR.stats'));
    values = [left; right];
    names = [left_names; right_names];
    index = zeros(87, 1);
    for roi = 1:87
        match = find(strcmp(names, roi_names{roi}));
        assert(numel(match) == 1, 'Each HFR ROI must occur exactly once across hemispheres.');
        index(roi) = match;
    end
    stats_data = struct();
    for f = 1:numel(fields)
        stats_data.(fields{f}) = values(index, f);
    end
    stats_data.region_names = roi_names;
    subject_out = fullfile(out_root, sub);
    if ~exist(subject_out, 'dir'), mkdir(subject_out); end
    save(fullfile(subject_out, 'morphology_stats_data.mat'), 'stats_data', '-v7');
end
