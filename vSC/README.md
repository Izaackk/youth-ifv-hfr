# vSC: structural-connectome variability

This module contains the selected code for computing ROI-wise inter-subject variability of individualized HFR structural-connectome matrices.

## Definition

For each age-sorted sliding window, the connectivity profile of each ROI is compared between every pair of subjects. Structural-connectome variability is defined as:

```text
vSC = 1 - mean pairwise connectivity-profile correlation
```

Analysis settings:

- Input connectome: individualized HFR88 structural-connectome matrix.
- ROI handling: remove the first non-cortical label, leaving 87 cortical ROIs.
- Window size: 20 subjects.
- Step size: 10 subjects.
- Window covariates: mean age, sex, and mean framewise displacement.

## Code

- `code/01_structural_connectome/SC_native2HFRmatrix_Parallel_withlog.sh`
  - Converts individualized surface annotations to DWI-space HFR88 parcellations.
  - Uses an existing SIFT-filtered tractogram to generate one structural-connectome matrix per subject.
  - Accepts a subject-list file as the first argument and the number of parallel jobs as the second argument.
- `code/02_sliding_window/build_vsc_win20_step10.m`
  - Loads the individual HFR88 matrices and computes win20/step10 ROI-wise vSC.
  - Writes `data_bin20_step10.mat` containing `var_bins`, `age_bins`, `sex_bins`, `mFD_bins`, and `wb_mean`.
- `code/utils/my_matcorr.m`
  - Computes matched column-wise correlations used for ROI connectivity-profile similarity.

## Usage

Replace paths beginning with `.../` before running the scripts.

Build individual structural-connectome matrices:

```bash
bash code/01_structural_connectome/SC_native2HFRmatrix_Parallel_withlog.sh \
  .../subject_ids_without_V1_MR.txt 10
```

Then configure the input paths at the top of `build_vsc_win20_step10.m` and run it in MATLAB.

## Requirements

- FreeSurfer
- MRtrix3
- ANTs
- GNU parallel
- MATLAB
- Individualized HFR annotation files
- HFR label-conversion LUT files
- Precomputed DWI-to-T1 affine transforms, mean b0 images, and `sift_1M.tck` tractograms

## Path and subject conventions

- The shell subject list contains IDs without `_V1_MR`, for example `HCD0001305`.
- The MATLAB subject list contains IDs with `_V1_MR`, for example `HCD0001305_V1_MR`.
- Individual matrices are stored as `<ID>_HFR_88_indi_SC.csv` or `.mat` inside a directory named `<ID>`.

## Verification

The archived win20/step10 input contained 57 windows and 87 ROIs. Formula-level checks reproduced the archived vSC significance and vSC-vIFP overlap masks exactly.
