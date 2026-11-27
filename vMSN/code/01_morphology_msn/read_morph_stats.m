function [values, names] = read_morph_stats(filename)
% Accept numeric-first mris_anatomical_stats stdout and name-first table rows.
lines = splitlines(string(fileread(filename)));
values = zeros(0, 9);
names = cell(0, 1);
for i = 1:numel(lines)
    line = strtrim(lines(i));
    if strlength(line) == 0 || startsWith(line, '#'), continue; end
    tokens = textscan(char(line), '%s');
    parts = string(tokens{1});
    if numel(parts) ~= 10, continue; end
    if isfinite(str2double(parts(1)))
        row = str2double(parts(1:9));
        name = parts(10);
    else
        row = str2double(parts(2:10));
        name = parts(1);
    end
    if ~all(isfinite(row)), continue; end
    values(end + 1, :) = row'; %#ok<AGROW>
    names{end + 1, 1} = char(name); %#ok<AGROW>
end
assert(~isempty(values), 'No morphometric data rows found in %s.', filename);
end
