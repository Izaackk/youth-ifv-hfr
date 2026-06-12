# vSC target verification

Target folder:

`.../HFR_output/MatchRate1/SC_DWI_20260102/GAM/k3`

Input MAT:

`.../HFR_output/MatchRate1/SC_DWI_20260102/Inter_new/data_bin20_step10.mat`

Observed variables:

- `var_bins`: 57 windows x 87 ROIs
- `age_bins`: 57 x 1
- `sex_bins`: 57 x 1
- `mFD_bins`: 57 x 1

Post-GAM verification:

- FDR source column: `GAM_Result.csv::anova_pvalue`
- BH-FDR threshold: `< 0.05`
- vSC significant ROIs: 57
- vIFP significant ROIs: 54
- vSC-vIFP overlap ROIs: 36
- Exact match to archived `deltaR2map_pFDR0.05/ROI_sig.mat`: true
- Exact match to archived `deltaR2map_SCagenocovmFD_cmpdiceage_yrk3GAMafdrp0.05/logi_ROI_selected.mat`: true
