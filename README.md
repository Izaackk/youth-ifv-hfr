# YouthIFV-HFR

Code accompanying the manuscript:

**Spatiotemporal Pattern and Structural Substrates of Individual Functional Variability in Youth**

This repository provides selected analysis code used to quantify inter-individual variability in youth brain organization and to model its age-related patterns. The current public subset focuses on individualized functional parcellation variability, structural-connectome variability, functional-connectivity variability, and generalized additive models.

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
```

## Analysis Modules

### vIFP

`vIFP/` contains code for computing variability of individualized functional parcellation. Pairwise ROI overlap is quantified with Dice overlap, and vIFP is defined as:

```text
vIFP = 1 - Dice overlap
```

The module includes scripts to compute pairwise cross-subject matched-ROI Dice overlap and to build age-sorted sliding-window vIFP summaries.

### vSC

`vSC/` contains code for individualized HFR structural-connectome construction and sliding-window structural-connectome variability. For each window, ROI connectivity profiles are compared across subject pairs, and vSC is defined as:

```text
vSC = 1 - mean pairwise connectivity-profile correlation
```

The structural-connectome workflow assumes precomputed diffusion preprocessing, DWI-to-T1 transforms, and SIFT-filtered tractograms.

### vFC

`vFC/` contains scripts for within-subject and between-subject functional-connectivity variability based on HFR ROI-to-ROI functional-connectivity matrices. The scripts compute intra-subject variability across sessions and inter-subject variability within age-sorted sliding windows.

### GAMs

`GAMs/` contains R scripts for generalized additive models of age-related effects. The scripts are written as reusable templates for ROI-wise or whole-brain average measures.

The main model form is:

```text
Measure ~ s(Age, bs = "cs", k = 3) + Sex + mFD
```

with a covariate-only null model used to estimate age-related delta R-squared.

## Typical Workflow

1. Compute pairwise vIFP from individualized matched ROI masks.
2. Build age-sorted sliding-window vIFP summaries.
3. Build individualized HFR structural-connectome matrices from tractography outputs.
4. Compute sliding-window vSC from HFR structural-connectome matrices.
5. Compute intra- and inter-subject vFC from HFR functional-connectivity matrices.
6. Run GAM scripts on the resulting window-level or subject-level matrices.

## Data Requirements

This repository contains code only. Input data, intermediate matrices, and manuscript result outputs are not included.

Expected inputs include:

- HCP-D subject lists and demographic/covariate tables.
- Individualized HFR ROI/parcellation files.
- Pairwise or subject-level functional-connectivity matrices.
- Individualized HFR structural-connectome matrices or the diffusion inputs needed to build them.
- Surface or ROI-level measures used as GAM inputs.

Because some source datasets are controlled-access or too large for GitHub, users should replace the placeholder paths with their own local data locations.

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

- MATLAB
- R
- R packages: `R.matlab`, `mgcv`, `gratia`, `visreg`, `ggplot2`, `dplyr`, `tidyr`
- FreeSurfer
- MRtrix3
- ANTs
- GNU parallel

Some scripts also call local helper functions included in the relevant module directories.

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
