% function intra_var(DatPAath, OutPath, subs)

clear,clc;

DataPath = '.../MatchRate1/FC_HFR_multi_sessions';
PathOut = '.../MatchRate1/Variability_FC';

subs = textread(['.../sublist.txt'], '%s');

list_sess = {'surf1_AP', 'surf1_PA', 'surf2_AP', 'surf2_PA'};
list_template = {'Atlas', 'Indi'};

num_HFR = 87;
IntraVariance = zeros(numel(subs), num_HFR);

% logi_nets = ones(1, num_HFR);
% logi_nets(:, 39: 42) = 0;
% logi_nets = logical(logi_nets);

for ntem = 1: 2
    tem = list_template{ntem};
    OutPath = fullfile(PathOut, tem, 'Intra_org');
    mkdir(OutPath)
    for s = 1: numel(subs)
        sub = subs{s}
        for nsess = 1: numel(list_sess)
            sess = list_sess{nsess};
            load(fullfile(DataPath, ['FC_', sess], ['ROI2ROIFC_', tem], [sub, '_big_corr.mat']))
%             FC(isnan(FC)) = 0;
%             Rmat(nsess, :, :) = FC(logi_nets, logi_nets);
            Rmat(nsess, :, :) = CorrMat;
        end

        Rmat(isnan(Rmat)) = 0;
%         count = 0;
        AveRmat = zeros(num_HFR, 1);

        for nroi = 1: num_HFR
            Rmat_roi = squeeze(Rmat(:, :, nroi));
            Rmat_roi(:, nroi) = []; 
            sim_intra_sub = corr(Rmat_roi');
            AveRmat(nroi) = mean(sim_intra_sub(tril(true(size(sim_intra_sub)), -1))); 
        end

%         IntraVariance(s, logi_nets) = 1 - AveRmat;
        IntraVariance(s, :) = 1 - AveRmat;

%         for m = 1: numel(list_sess)
%             for n = m + 1: numel(list_sess)
%                 count = count + 1;
%                 tmp = my_matcorr(squeeze(Rmat(m, :, :))', squeeze(Rmat(n, :, :))');
%                 tmp(isnan(tmp)) = 0;
%                 %             AveRmat = AveRmat + diag(tmp);
%                 AveRmat = AveRmat + tmp';
%             end
%         end
%         count
%         IntraVariance(s, :) = 1 - AveRmat./count;

        path_save = [OutPath, '/', sub, '/'];
        mkdir(path_save)

        save_mgh_HFR_fs5(IntraVariance(s, :), [path_save '/intravariance_HFR87'])

        intravariance = IntraVariance(s, :);
        save([path_save '/intravariance_HFR87/intravariance_HFR87.mat'], 'intravariance')
    end
    meanIntraVariance = mean(IntraVariance);
    save_mgh_HFR_fs5(meanIntraVariance, [OutPath '/meanIntravariance_across4sess_HFR87'])

    save([OutPath '/meanIntravariance_across4sess_HFR87/meanIntravariance_across4sess_HFR87.mat'], 'meanIntraVariance')
end