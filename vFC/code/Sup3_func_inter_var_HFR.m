function Sup3_func_inter_var_HFR(DataPath, OutPath, IntraPath, subs, template)

fprintf('\n=== Start Sup_func_inter_var_HFR87  | Template: %s ===\n', template);
fprintf('IntraPath: %s\n', IntraPath);

mkdir(OutPath);

num_HFR = 87;
list_sess = {'surf1_AP', 'surf1_PA', 'surf2_AP', 'surf2_PA'};

Sim_allsession = zeros(numel(list_sess), num_HFR);
Var_allsession = zeros(numel(list_sess), num_HFR);

%% 
for i = 1:numel(list_sess)
    sess = list_sess{i};
    fprintf('  -> Processing session %d / %d : %s\n', i, numel(list_sess), sess);

    Rmat = nan(numel(subs), num_HFR, num_HFR);
    for s = 1:numel(subs)
        sub = subs{s};
        f_FC = fullfile(DataPath, ['FC_', sess], ['ROI2ROIFC_', template], [sub, '_big_corr.mat']);
        if exist(f_FC, 'file')
            load(f_FC, 'CorrMat');
            Rmat(s, :, :) = CorrMat;
        else
            warning('Missing FC file: %s', f_FC);
        end
    end

    Rmat(isnan(Rmat)) = 0;

    AveRmat = zeros(num_HFR, 1);
    for nroi = 1:num_HFR
        Rmat_roi = squeeze(Rmat(:, :, nroi));
        Rmat_roi(:, nroi) = [];
        sim_inter = corr(Rmat_roi');
        AveRmat(nroi) = mean(sim_inter(tril(true(size(sim_inter)), -1)));
    end

    InterSimilarity_session = AveRmat;
    InterVariance_session   = 1 - AveRmat;

    path_save = fullfile(OutPath, ['session', num2str(i)]);
    mkdir(path_save);

    save_mgh_HFR_fs5(InterSimilarity_session, fullfile(path_save, 'similarity_HFR87'));
    save(fullfile(path_save, 'similarity_HFR87', 'InterSimilarity.mat'), 'InterSimilarity_session');

    save_mgh_HFR_fs5(InterVariance_session, fullfile(path_save, 'intervariance_HFR87'));
    save(fullfile(path_save, 'intervariance_HFR87', 'InterVariance.mat'), 'InterVariance_session');

    Sim_allsession(i, :) = InterSimilarity_session;
    Var_allsession(i, :) = InterVariance_session;
end

%% 
meanInterSimilarity = mean(Sim_allsession, 1);
meanInterVariance   = mean(Var_allsession, 1);

save_mgh_HFR_fs5(meanInterSimilarity, fullfile(OutPath, 'meanSimilarity_across4sess_HFR87'));
save(fullfile(OutPath, 'meanSimilarity_across4sess_HFR87', 'meanInterSimilarity.mat'), 'meanInterSimilarity');

save_mgh_HFR_fs5(meanInterVariance, fullfile(OutPath, 'meanInterVariance_across4sess_HFR87'));
save(fullfile(OutPath, 'meanInterVariance_across4sess_HFR87', 'meanInterVariance.mat'), 'meanInterVariance');

%% 
f_intra_win = fullfile(OutPath, 'meanIntraVariance_win.mat');
if exist(f_intra_win, 'file')
    load(f_intra_win, 'meanIntra_win');
    Intra_value = meanIntra_win(:);
else
    warning('meanIntraVariance_win.mat not found at %s, skipping regression step.', f_intra_win);
    Intra_value = zeros(num_HFR, 1);
end

Variability = zeros(numel(list_sess), num_HFR);        
Variability_norm = zeros(numel(list_sess), num_HFR);  

for i = 1:numel(list_sess)
    Inter = Var_allsession(i, :)';
    % Inter = β₁·Intra + β₀ + ε
    X = [Intra_value, ones(num_HFR, 1)];
    beta = pinv(X) * Inter;

    % Variability_norm
    tmp_resid = Inter - X * beta;
    Variability_norm(i, :) = tmp_resid';

    % Variability
    tmp_keepIntercept = Inter - X(:,1) * beta(1);
    Variability(i, :) = tmp_keepIntercept';

    % 
    path_save = fullfile(OutPath, ['session', num2str(i)]);

    % 
    save_mgh_HFR_fs5(Variability_norm(i, :), fullfile(path_save, 'intervariability_norm_HFR87'));
    save(fullfile(path_save, 'intervariability_norm_HFR87', 'Variability_norm.mat'), 'Variability_norm');

    % 
    save_mgh_HFR_fs5(Variability(i, :), fullfile(path_save, 'intervariability_HFR87'));
    save(fullfile(path_save, 'intervariability_HFR87', 'Variability.mat'), 'Variability');
end

% 
meanVariability_norm = mean(Variability_norm, 1);
meanVariability = mean(Variability, 1);

% 
save_mgh_HFR_fs5(meanVariability_norm, fullfile(OutPath, 'InterVariability_norm_HFR87'));
save(fullfile(OutPath, 'InterVariability_norm_HFR87', 'meanVariability_norm.mat'), 'meanVariability_norm', 'Variability_norm');

save_mgh_HFR_fs5(meanVariability, fullfile(OutPath, 'InterVariability_HFR87'));
save(fullfile(OutPath, 'InterVariability_HFR87', 'meanVariability.mat'), 'meanVariability', 'Variability');

fprintf('=== Sup_func_inter_var_HFR87 finished successfully ===\n');
% end
