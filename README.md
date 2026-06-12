# YouthIFV-HFR

Code accompanying the manuscript:

**Spatiotemporal Pattern and Structural Substrates of Individual Functional Variability in Youth**

This repository currently contains the selected code for:

- `vIFP/`: computing inter-subject variability of individualized functional parcellation, defined as `1 - Dice overlap`, and building sliding-window vIFP summaries.
- `vSC/`: generating individualized HFR structural-connectome matrices and computing win20/step10 structural-connectome variability.
- `GAMs/`: generalized additive model scripts for age-related effects and whole-brain average curves.

## Path Placeholders

Paths beginning with `.../` are placeholders. Replace them with the corresponding local data or output root before running the scripts.

## Data

This repository contains code only. Input data, intermediate matrices, and manuscript result outputs are not included.

The structural-connectome workflow additionally requires FreeSurfer, MRtrix3, ANTs, GNU parallel, individualized HFR annotations, and precomputed tractography inputs.

## Notes

The released code is organized to expose the computational method and main analysis workflow. Historical exploratory or sensitivity scripts are not included in this initial public subset.
