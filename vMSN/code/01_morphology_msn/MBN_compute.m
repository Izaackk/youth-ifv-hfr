clear; clc;

% Build Pearson MSN matrices from eight min-max-normalized ROI features.
data_path = '.../HFR_output/MatchRate1/MBC/data_HFR';
sid_file = '.../SubList.txt';
out_root = '.../HFR_output/MatchRate1/MBC/data_HFR';

SIDs = strtrim(string(splitlines(fileread(sid_file))));
SIDs(SIDs == "") = [];
feature_names = {'surface_area', 'gray_matter_vol', 'thickness_mean', ...
    'thickness_std', 'mean_curv', 'gauss_curv', 'folding_index', 'curv_index'};
roi_names = cellstr(compose('Network_%d', (2:88)'));
for n = 1:numel(SIDs)
    sub = char(SIDs(n));
    filename = fullfile(data_path, sub, 'morphology_stats_data.mat');
    if ~isfile(filename)
        % Compatibility with the original filename containing an embedded space.
        filename = fullfile(data_path, sub, 'morphology_ stats_data.mat');
    end
    loaded = load(filename, 'stats_data');
    stats_data = loaded.stats_data;
    index = zeros(87, 1);
    for roi = 1:87
        match = find(strcmp(stats_data.region_names, roi_names{roi}));
        assert(numel(match) == 1, 'Each selected HFR ROI must have exactly one stats row.');
        index(roi) = match;
    end
    mph_feature = zeros(8, 87);
    for f = 1:8
        values = stats_data.(feature_names{f});
        mph_feature(f, :) = reshape(values(index), 1, []);
    end
    assert(all(isfinite(mph_feature(:))), 'Morphometric features must be finite.');
    % Normalize each feature across ROIs, not each ROI across features.
    range_values = max(mph_feature, [], 2) - min(mph_feature, [], 2);
    range_values(range_values == 0) = 1;
    mph_feature = (mph_feature - min(mph_feature, [], 2)) ./ range_values;
    MBC = corr(mph_feature);
    subject_out = fullfile(out_root, sub);
    if ~exist(subject_out, 'dir'), mkdir(subject_out); end
    save(fullfile(subject_out, 'MBC_minmaxnorm_novertices.mat'), ...
        'MBC', 'roi_names', 'feature_names', '-v7');
end
