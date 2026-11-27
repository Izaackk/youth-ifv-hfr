# YouthIFV-HFR

Code accompanying the manuscript:

**Spatiotemporal Pattern and Structural Substrates of Individual Functional Variability in Youth**

This repository provides selected analysis code used to quantify inter-individual variability in youth brain organization and to model its age-related patterns. The current public subset includes individualized functional parcellation, structural-connectome, functional-connectivity, morphometric-similarity-network, and myelin-distribution variability, together with generalized additive models.

## Repository Contents

```text
GAMs/
  GAM_izaac_version2.R
  GAM_izaac_multiWB.R

vIFP/
  code/01_pairwise_cmpdice/
  code/02_sliding_window_vifp/

vSC/
  code/01_structural_connectome/
  code/02_sliding_window/
  code/utils/

vFC/
  code/

vMSN/
  code/01_morphology_msn/
  code/02_sliding_window/

vMD/
  code/01_myelin_roi/

utils/
  Fisherz.m
  matcorr.m
  save_mgh_HFR_fs5.m
```

## Analysis Modules

### vIFP

`vIFP/` contains code for computing variability of individualized functional parcellation. Pairwise ROI overlap is quantified with Dice overlap, and vIFP is defined as:

```text
vIFP = 1 - Dice overlap
```

The module includes scripts to compute pairwise cross-subject matched-ROI Dice overlap and to build age-sorted sliding-window vIFP summaries.

Run `dice_HFR_group.m` first. Its `vIFP/outputs/CmpDice/cmp_dice_results_new.mat` contains **vIFP, not Dice overlap**, together with subject IDs and an explicit measure marker. Point `pairwise_file` in the subgroup and sliding-window scripts to this file. They use its subject order and must not apply `1 - x` again. Unmarked historical files with the same filename can contain Dice overlap and are rejected; regenerate the marked input from matched ROI masks. Pairwise means omit undefined values when both ROI masks are empty.

The separate `All_601/dice.mat` retains legacy Dice overlap; `All_601/vifp.mat` contains its complementary vIFP measure. The subject list supplied to the pairwise builder must match the order of `All_ROI_Ind_Big`.

### vSC

`vSC/` contains code for individualized HFR structural-connectome construction and sliding-window structural-connectome variability. For each window, ROI connectivity profiles are compared across subject pairs, and vSC is defined as:

```text
vSC = 1 - mean pairwise connectivity-profile correlation
```

The structural-connectome workflow assumes precomputed diffusion preprocessing, DWI-to-T1 transforms, and SIFT-filtered tractograms.

Edit paths, window settings, and optional excluded IDs at the top of `build_vsc_win20_step10.m`. It reads 88 x 88 MAT or CSV matrices without writing to the input directory, zeroes the diagonal, and removes the first label before profile comparison. Zero-variance profile correlations are replaced with zero, matching the original calculation.

### vFC

`vFC/` contains scripts for within-subject and between-subject functional-connectivity variability based on HFR ROI-to-ROI functional-connectivity matrices. The scripts compute intra-subject variability across sessions and inter-subject variability within age-sorted sliding windows.

Run the intra-subject script before the inter-subject script. Surface exports use the included `utils/save_mgh_HFR_fs5.m`; configure its `template_dir` and supply the authorized HFR87-to-116 mapping and fsaverage5 label surfaces. These templates are not distributed. An intra-subject age model uses one row per subject, not sliding-window averages; assembling its covariate table is a separate input-preparation step.

### vMSN

`vMSN/` extracts eight morphometric features from native HFR statistics, min-max normalizes each feature across 87 ROIs, and computes subject-level Pearson morphometric-similarity networks. Window-level variability is `1 - mean pairwise correlation` of Fisher-z-transformed connectivity profiles. The self-connection is retained, matching the original analysis.

### vMD

`vMD/` computes the Kruskal-Wallis H statistic across subjects' vertex-level myelin distributions within each individualized ROI. Its whole-brain summary is the unweighted mean of the 87 ROI H values. This is not `1 - correlation` or a variance of ROI means.

Both modules default to age-sorted windows of 20 subjects with a step of 10. Supply the subject list for the intended analysis; the original 601-subject references have 59 complete windows. A 569-subject multimodal list gives 55 windows and is a different analysis.

Their main analysis files are scripts with editable paths and parameters at the top. Supporting functions handle file parsing and metric calculations. Configure and run the scripts directly with the inputs and software listed below.

### GAMs

`GAMs/` contains R scripts for generalized additive models of age-related effects. The scripts are written as reusable templates for ROI-wise or whole-brain average measures.

The main model form is:

```text
Measure ~ s(Age, bs = "cs", k = 3) + Sex + mFD
```

with a covariate-only null model used to estimate age-related delta R-squared. This describes `GAM_izaac_version2.R`. The current `GAM_izaac_multiWB.R` fits age plus sex without mFD in its GAMs; it is not an interchangeable covmFD reproducer. Select the intended model formula before fitting.

The two GAM scripts retain their previously published implementation. The whole-brain script fits the single-column summary selected in each `datasets` entry; it does not average an ROI matrix. Configure both the input file and variable names for the chosen module. `R.matlab` converts underscores in MATLAB variable names to dots:

| Module | ROI measure | Sex covariate | Whole-brain measure |
| --- | --- | --- | --- |
| vIFP | `vifp.bins` | `sex.de.bins` | `wb.mean.vifp` in the same window file |
| vSC | `var.bins` | `sex.bins` | `wb.mean` in the same window file |
| vFC inter-subject | `var.bins` | `sex.bins` | Prepare a single-column ROI mean separately |
| vMSN | `H` or `var.all` | `sex.de.bins` or `sex.bins` | `mvar.all` in the same window file |
| vMD | `H` | `sex.de.bins` | `H` in the separate `WB` file |

All window inputs also require `age.bins` and `mFD.bins`. Edit the ROI script's input assignments as needed; its defaults target vMD. The ROI script also needs a user-supplied SA-axis file with one value per ROI. Its derivative grid is fixed to ages 8-22 and its plot limits are measure-specific; review both for other datasets.

`R_sq` is adjusted R-squared; `delta_R_sq` is the full-minus-null adjusted R-squared difference with a sign assigned from the mean age derivative. `p_value` is the smooth-term p-value, whereas `anova_pvalue` comes from the legacy two-model comparison. The `anova_p0.05` flag is uncorrected, not FDR-adjusted. The shrinkage smooth (`bs = "cs"`) makes that model-comparison test potentially anti-conservative; see the [mgcv ANOVA documentation](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/anova.gam.html). Overlapping windows are not independent observations. These scripts preserve the original methods, not a newly validated inference procedure; neither reverting the scripts nor a successful software test resolves these statistical limitations.

## Typical Workflow

1. Compute pairwise vIFP from individualized matched ROI masks.
2. Build age-sorted sliding-window vIFP summaries.
3. Build individualized HFR structural-connectome matrices from tractography outputs.
4. Compute sliding-window vSC from HFR structural-connectome matrices.
5. Compute intra- and inter-subject vFC from HFR functional-connectivity matrices.
6. Extract morphometric features, build MSN matrices, and compute sliding-window vMSN.
7. Compute sliding-window vMD from aligned myelin surfaces and individualized parcellations.
8. Run GAM scripts on the resulting window-level or subject-level matrices.

## Data Requirements

This repository contains code only. Input data, intermediate matrices, and manuscript result outputs are not included.

Expected inputs include:

- HCP-D subject lists and demographic/covariate tables.
- Individualized HFR ROI/parcellation files.
- Pairwise or subject-level functional-connectivity matrices.
- Individualized HFR structural-connectome matrices or the diffusion inputs needed to build them.
- Surface or ROI-level measures used as GAM inputs.
- Native HFR morphometric statistics for vMSN.
- fsaverage4 myelin maps and a cortical vertex-index mask for vMD.
- HFR label-conversion LUTs for structural-connectome construction and the matched-ROI mapping/label surfaces for vFC surface exports.

Because some source datasets are controlled-access or too large for GitHub, users should replace the placeholder paths with their own local data locations.

Access to research data does not by itself authorize redistribution. Obtain controlled inputs through their authorized source. The new vMSN/vMD modules contain no participant data or private template files.

## Path Placeholders

Paths beginning with `.../` are placeholders. Replace them before running the scripts, for example:

```matlab
ROIPath = '.../HFR_output/MatchRate1';
SubIDs = textread('.../SubList.txt', '%s');
```

The placeholders are intentionally kept visible so that users can adapt the scripts to their own data organization.

## Software Requirements

The full workflow uses a mixture of MATLAB, R, and neuroimaging command-line tools.

Core requirements:

- MATLAB (tested with R2022b) and Statistics and Machine Learning Toolbox
- R (tested with 4.4.1)
- R packages: `R.matlab`, `mgcv`, `gratia`, `visreg`, `ggplot2`, `paletteer`, `scales`, `dplyr`, `pals`, `tidyr`
- FreeSurfer
- MRtrix3
- ANTs
- GNU parallel

The repository-root `utils/` directory includes MATLAB helpers used by vMSN and vFC:

- `utils/Fisherz.m`: the finite Fisher transformation used for MSN entries.
- `utils/matcorr.m`: correlation between corresponding matrix columns.
- `utils/save_mgh_HFR_fs5.m`: original matched-ROI mapping for vFC surface exports; requires user-supplied templates and FreeSurfer `load_mgh`/`save_mgh`.

The vMSN and vMD main scripts read subject information directly and do not depend on `shared/`. Supply a text subject list and a covariate table with `SubID`, numeric `Age`, numeric `Sex`, and numeric `meanFD`; rows are matched to the subject-list order. Local tests and module-level README files are excluded from publication.

vMD additionally needs MATLAB readers `gifti` and FreeSurfer `load_mgh` on the MATLAB path. Input surfaces must already be aligned; these modules do not perform native parcellation generation or surface resampling.

The GAM environment used for local checks includes `mgcv` 1.9-1, `gratia` 0.9.2, and `visreg` 2.7.0. In that environment `gratia` warns that `term` in `derivatives()` is deprecated; the restored scripts have not been migrated to a different API.

## Verification Scope

- vSC: all saved fields agree exactly before and after script cleanup for 582 subjects; the 57 x 87 window matrix and covariates match the archived reference exactly.
- vMSN: 601 archived morphometric feature files reproduce the subject MSN matrices, the 59 x 87 variability matrix, and window covariates exactly. Current native statistics differ from the archived feature inputs for 596 subjects, so raw-statistics-to-manuscript reproduction is not established.
- vMD: the 601-subject surface-input calculation reproduces the 59 x 87 ROI H matrix and window covariates exactly.
- Local synthetic tests cover vIFP measure propagation, vFC calculations and surface export, and MATLAB/R input-output workflows. Tests and their generated inputs remain local-only.

These checks do not establish a new reproduction of all manuscript GAM results or validate upstream preprocessing. No study data or private template is included in the repository.

## Outputs

Typical outputs include:

- Pairwise Dice/vIFP matrices.
- Sliding-window matrices such as `data_bin20_step10.mat`.
- ROI-wise variability matrices.
- GAM summary tables such as `GAM_Result.csv`.
- Fitted age curves and derivative summaries.

Large generated outputs are intentionally ignored by `.gitignore`.

## Notes

This repository is organized to expose the computational method and main analysis workflow. Historical exploratory scripts, intermediate result files, and large neuroimaging outputs are not included in this public subset.

The code is being prepared alongside the manuscript and may be expanded as additional analysis modules are finalized for release.
