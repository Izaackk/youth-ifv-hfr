% function InterSub_variance_slidingwindow_HFR87(window_size, step_size)

clear; clc;

%%
DataPath  = '.../MatchRate1/FC_HFR_multi_sessions';
BasePath  = '.../MatchRate1/Variability_FC';
InfoPath  = '.../info';

SIDs      = textread(fullfile(InfoPath, 'sublist.txt'), '%s');

Data_info = readtable(fullfile(InfoPath, 'Subinfo.xlsx'));

Allsubs   = Data_info.SubID;
Allage    = Data_info.Age;
Allsex    = Data_info.Sex;
AllmeanFD = Data_info.meanFD;

list_template = {'Atlas', 'Indi'};

%%
for ntem = 1: numel(list_template)
    tem = list_template{ntem};
    fprintf('\n===== Processing Template: %s =====\n', tem);

    IntraPath = fullfile(BasePath, tem, 'Intra_org');
    OutPathBase = fullfile(BasePath, tem, 'Inter');

    num_HFR = 87;

    %% 
    for nsub = 1: numel(SIDs)
        sub = SIDs{nsub};
        ind = ismember(Allsubs, sub);
        age(nsub, 1)    = Allage(ind);
        sex(nsub, 1)    = Allsex(ind);
        meanFD(nsub, 1) = AllmeanFD(ind);
    end    

    %% 
    [~, I] = sort(age);
    window_size = 20;
    step_size = 10;

    num_windows = floor((numel(age) - window_size) / step_size) + 1;
    var_bins = zeros(num_windows, num_HFR);
    age_bins = zeros(num_windows, 1);
    sex_bins = zeros(num_windows, 1);
    mFD_bins = zeros(num_windows, 1);

    %% 
    for nwin = 1: num_windows
        disp(['Window ', num2str(nwin), ' / ', num2str(num_windows)]);

        start_idx = (nwin - 1) * step_size + 1;
        end_idx = start_idx + window_size - 1;
        I_win = I(start_idx:end_idx);

        sub_win = SIDs(I_win);
        age_win = age(I_win);
        sex_win = sex(I_win);
        mFD_win = meanFD(I_win);

        age_bins(nwin, 1) = mean(age_win);
        sex_bins(nwin, 1) = mean(sex_win);
        mFD_bins(nwin, 1) = mean(mFD_win);

        % 
        OutPath = fullfile(OutPathBase, ...
            ['win', num2str(window_size), '_step', num2str(step_size)], ...
            ['win', num2str(nwin)]);
        mkdir(OutPath);

        % 
        intra_win = nan(numel(sub_win), num_HFR);
        for ns = 1:numel(sub_win)
            sub = sub_win{ns};
            f_intra = fullfile(IntraPath, sub, 'intravariance_HFR87', 'intravariance_HFR87.mat');
            if exist(f_intra, 'file')
                load(f_intra, 'intravariance');
                intra_win(ns, :) = intravariance;
            else
                warning('Missing intra file: %s', f_intra);
            end
        end
        meanIntra_win = mean(intra_win, 1);

        % 
        save(fullfile(OutPath, 'meanIntraVariance_win.mat'), 'meanIntra_win', 'sub_win', 'age_win', 'sex_win', 'mFD_win');

        % 
        Sup3_func_inter_var_HFR(DataPath, OutPath, IntraPath, sub_win, tem);

        % 
        load(fullfile(OutPath, 'InterVariability_HFR87', 'meanVariability.mat'));
        var_bins(nwin, :) = meanVariability;
    end

    %% 
    [r, p] = corr(var_bins, age_bins);

    num_permutations = 10000;
    count = zeros(size(r));
    r_perm_dstr = zeros(size(r, 1), num_permutations);

    for i = 1:num_permutations
        perm_age = age_bins(randperm(length(age_bins)));
        [r_perm, ~] = corr(var_bins, perm_age);
        r_perm_dstr(:, i) = r_perm;
        count(abs(r_perm) >= abs(r)) = count(abs(r_perm) >= abs(r)) + 1;
    end

    p_perm = count / num_permutations;

    save(fullfile(OutPathBase, ...
        ['win', num2str(window_size), '_step', num2str(step_size)], ...
        'corr_results_HFR87.mat'), ...
        'r', 'p', 'p_perm', 'r_perm_dstr', ...
        'var_bins', 'age_bins', 'sex_bins', 'mFD_bins');
end
