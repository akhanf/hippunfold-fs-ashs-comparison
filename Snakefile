import pandas as pd

df = pd.read_table('resources/subjects_HCA.csv',sep=',')


subjects, = glob_wildcards('freesurfer/sub-{subject}/mri/lh.hippoAmygLabels-T1.v21.CA.mgz')

rule all:
    input:
        expand('results/{method}_volumes.tsv',
            method=['freesurfer','ashs','hippunfold-magdeburgatlas','hippunfold-freesurferatlas','hippunfold-bigbrainatlas'])



rule mgz_to_nii:
    input:
        'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T1.v21.CA.mgz'
    output:
        'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T1.v21.CA.nii.gz'
    shell: 
        'mri_convert {input} {output}'


rule get_subfield_vols_fs:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T1.v21.CA.nii.gz',
                hemi=['lh','rh'], allow_missing=True),
        lookup_tsv='resources/freesurfer_v21_CA_subfields_dseg.tsv'
    output:
        tsv='results/freesurfer/sub-{subject}_volumes.tsv'
    script:
        "scripts/gen_volume_tsv.py"

rule get_subfield_vols_ashs:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'ashs/sub-{subject}/final/{subject}_{hemi}_lfseg_corr_nogray.nii.gz',
                hemi=['left','right'], allow_missing=True),
        lookup_tsv='resources/desc-subfields_atlas-magdeburg_dseg.tsv'
    output:
        tsv='results/ashs/sub-{subject}_volumes.tsv'
    script:
        "scripts/gen_volume_tsv.py"


rule get_subfield_vols_hippunfold:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'hippunfold_atlasfix/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-cropT1w_desc-subfields_atlas-{atlas}_dseg.nii.gz',
                hemi=['L','R'], allow_missing=True),
        lookup_tsv='resources/hippunfold_desc-subfields_atlas-{atlas}_dseg.tsv'
    output:
        tsv='results/hippunfold-{atlas}atlas/sub-{subject}_volumes.tsv'
    script:
        "scripts/gen_volume_tsv.py"





rule concat_subj_vols_tsv:
    """Concatenate all subject tsv files into a single tsv"""
    input:
        tsv=expand('results/{method}/sub-{subject}_volumes.tsv',
            subject=subjects,allow_missing=True)
    output:
        tsv='results/{method}_volumes.tsv'
    run:
        pd.concat([pd.read_table(in_tsv) for in_tsv in input]).to_csv(
            output.tsv, sep="\t", index=False
        )
