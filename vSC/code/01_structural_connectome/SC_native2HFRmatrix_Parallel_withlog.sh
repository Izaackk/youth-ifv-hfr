#!/usr/bin/env bash

set -euo pipefail

path_dwi=".../HCPD_DWI/MRtrix3"
path_out=".../HCPd_SC_HFR88indi"
path_annot=".../HFR_output/MatchRate1/Indi_HFR_parcellation"
path_lut=".../HFR_LUT"
path_recon=".../HCPD_freesurfer"

subject_list="${1:-.../HCPd_info/subject_ids_without_V1_MR.txt}"
jobs="${2:-10}"
atlas="HFR_88_indi"

log_dir="${path_out}/logs"
mkdir -p "${log_dir}"

process_subject() {
    local subj="$1"
    local fs_subject="${subj}_V1_MR"
    local subject_out="${path_out}/${fs_subject}"
    local log_file="${log_dir}/${subj}.log.txt"

    exec > >(tee -a "${log_file}") 2>&1

    echo "======================================"
    echo "Subject: ${subj}"
    echo "Start time: $(date)"
    echo "======================================"

    mkdir -p "${subject_out}"
    export SUBJECTS_DIR="${path_recon}"

    echo "[INFO] Copy individualized HFR annotations"
    cp "${path_annot}/${fs_subject}/Clusters/lh_HFR_88_native.annot" \
       "${path_recon}/${fs_subject}/label/lh.HFR_88_native.annot"
    cp "${path_annot}/${fs_subject}/Clusters/rh_HFR_88_native.annot" \
       "${path_recon}/${fs_subject}/label/rh.HFR_88_native.annot"

    echo "[INFO] Convert surface annotations to a volume"
    mri_aparc2aseg \
      --s "${fs_subject}" \
      --annot HFR_88_native \
      --o "${subject_out}/${fs_subject}_${atlas}.mgz"

    echo "[INFO] Convert the parcellation to MRtrix format"
    mrconvert -datatype uint32 \
      "${subject_out}/${fs_subject}_${atlas}.mgz" \
      "${subject_out}/${fs_subject}_${atlas}_native.mif"

    echo "[INFO] Reindex labels to 1-88"
    labelconvert \
      "${subject_out}/${fs_subject}_${atlas}_native.mif" \
      "${path_lut}/LUT_orig.txt" \
      "${path_lut}/LUT_order.txt" \
      "${subject_out}/${fs_subject}_${atlas}_native.coreg.mif"

    echo "[INFO] Convert the reindexed parcellation to NIfTI"
    mrconvert -datatype uint32 \
      "${subject_out}/${fs_subject}_${atlas}_native.coreg.mif" \
      "${subject_out}/${fs_subject}_${atlas}_native.coreg.mgz"

    mri_convert \
      "${subject_out}/${fs_subject}_${atlas}_native.coreg.mgz" \
      "${subject_out}/${fs_subject}_${atlas}_native.coreg.nii.gz"

    echo "[INFO] Transform the parcellation to DWI space"
    antsApplyTransforms -d 3 \
      -i "${subject_out}/${fs_subject}_${atlas}_native.coreg.nii.gz" \
      -r "${path_dwi}/${subj}_temp/mean_b0.nii.gz" \
      -n NearestNeighbor \
      -t "[${path_dwi}/${subj}_temp/ants0GenericAffine.mat,1]" \
      -o "${subject_out}/${fs_subject}_${atlas}_native_inDWI.nii.gz"

    echo "[INFO] Build the individualized structural-connectome matrix"
    tck2connectome \
      -symmetric -zero_diagonal -scale_invnodevol \
      -assignment_radial_search 2 \
      "${path_dwi}/${subj}_temp/sift_1M.tck" \
      "${subject_out}/${fs_subject}_${atlas}_native_inDWI.nii.gz" \
      "${subject_out}/${fs_subject}_${atlas}_SC.csv" \
      -out_assignment "${subject_out}/${fs_subject}_${atlas}_SC_assignment.csv" \
      -force

    echo "End time: $(date)"
    echo "Finished subject: ${subj}"
}

export -f process_subject
export path_dwi path_out path_annot path_lut path_recon atlas log_dir

parallel -j "${jobs}" process_subject :::: "${subject_list}"
