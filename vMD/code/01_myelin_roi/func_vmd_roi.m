function H = func_vmd_roi(myelin_window, parc_window, n_roi)
% Kruskal-Wallis H for vertex-level myelin values grouped by subject.
% Inputs have dimensions retained vertices x subjects.
assert(isequal(size(myelin_window), size(parc_window)), ...
    'Myelin and parcellation matrices must have identical dimensions.');
assert(size(myelin_window, 2) >= 2, 'At least two subjects are required.');
H = zeros(1, n_roi);
for roi = 1:n_roi
    samples = cell(size(myelin_window, 2), 1);
    labels = cell(size(samples));
    for s = 1:numel(samples)
        selected = parc_window(:, s) == roi;
        samples{s} = myelin_window(selected, s);
        labels{s} = s .* ones(sum(selected), 1);
    end
    values = vertcat(samples{:});
    groups = vertcat(labels{:});
    assert(numel(unique(groups(~isnan(values)))) >= 2, ...
        'Each ROI needs nonmissing vertices from at least two subjects.');
    [~, table] = kruskalwallis(values, groups, 'off');
    H(roi) = table{2, 5};
end
end
