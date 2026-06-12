# vIFP definition and reference check

Target reference:

- `.../HFR_output/MatchRate1/CmpDice/All_601/dice.mat`

Observed reference structure:

- Variable: `averages`
- Size: `87 x 1`
- Mean: `0.2802065665491873`

Audit result:

- `averages` in the reference file equals the lower-triangle mean of `cmp_dice_results_new.mat` for each ROI.
- Therefore the reference file stores Dice overlap, not vIFP.
- Corrected vIFP/cmpdice is `1 - Dice overlap`.
- Formula-level check on 2026-05-27:
  - Dice-overlap reproduction max absolute difference vs reference: `2.775557561562891e-16`
  - Dice-overlap mean: `0.2802065665491874`
  - Corrected vIFP mean: `0.7197934334508125`
  - Corrected vIFP range: `0.5303690944276207` to `0.9023302958434632`

Implementation in this folder:

- `code/01_pairwise_cmpdice/dice_HFR_group.m` recomputes pairwise Dice overlap from matched ROI masks and writes both legacy Dice and corrected vIFP outputs.
- `code/02_sliding_window_vifp/dice_slidingwindow.m` computes sliding-window vIFP as `1 - Dice overlap`.
