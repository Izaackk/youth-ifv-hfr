# vIFP / cmpdice

This folder contains the selected vIFP/cmpdice code for the manuscript code release.

Scope:

- Reproduce the historical `CmpDice/All_601/dice.mat` ROI map from matched-ROI Dice overlap.
- Compute corrected inter-subject vIFP/cmpdice as `1 - Dice overlap`.
- Build sliding-window vIFP summaries from `cmp_dice_results_new.mat`.
- Provide the win20/step10 vIFP input used by the age-modeling scripts.

Original code was not edited. Files under `code/` are copied from `.../HFR`.

## Structure

- `code/01_pairwise_cmpdice/`
  - Pairwise cross-subject matched-ROI Dice computation and All_601 map generation.
- `code/02_sliding_window_vifp/`
  - Sliding-window vIFP construction from the pairwise cmpdice matrix.
- `docs/selected_vifp_files.txt`
  - Source-to-destination manifest for this step.

## Important definition

The existing reference file `.../HFR_output/MatchRate1/CmpDice/All_601/dice.mat` stores mean Dice overlap in the variable `averages`.

For this cleaned module:

- `dice.mat`: legacy reproduction of mean Dice overlap.
- `vifp.mat`: corrected vIFP/cmpdice, where `averages = 1 - Dice overlap`.
- Sliding-window scripts now use vIFP (`1 - Dice overlap`) for `vifp_bins`.
- `dice_bins` is kept as an alias of `vifp_bins` only for compatibility with existing R scripts.

## Path placeholders

Paths beginning with `.../` are placeholders. Replace them with the corresponding local data root before running the scripts.
