import pandas as pd

df = pd.read_table('resources/subjects_HCA.csv',sep=',')


subjects, = glob_wildcards('freesurfer/sub-{subject}/mri/lh.hippoAmygLabels-T2.v21.CA.mgz')

rule all:
    input:
        expand('results/{method}_volumes.tsv',
            method=['freesurfer','ashs',
                    'hippunfold-t1-magdeburgatlas','hippunfold-t1-freesurferatlas','hippunfold-t1-bigbrainatlas',
                    'hippunfold-t2-magdeburgatlas','hippunfold-t2-freesurferatlas','hippunfold-t2-bigbrainatlas'])



rule mgz_to_nii:
    input:
        'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.mgz'
    output:
        'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.nii.gz'
    shell: 
        'mri_convert {input} {output}'

rule all_fs_registered:
    input:
        segs=expand(
            'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.space-T2w.nii.gz',
                hemi=['lh','rh'], subject=subjects,allow_missing=True),
 

rule get_subfield_vols_fs:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.nii.gz',
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


rule get_subfield_vols_hippunfold_t1:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'hippunfold_atlasfix/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-cropT1w_desc-subfields_atlas-{atlas}_dseg.nii.gz',
                hemi=['L','R'], allow_missing=True),
        lookup_tsv='resources/hippunfold_desc-subfields_atlas-{atlas}_dseg.tsv'
    output:
        tsv='results/hippunfold-t1-{atlas}atlas/sub-{subject}_volumes.tsv'
    script:
        "scripts/gen_volume_tsv.py"

rule get_subfield_vols_hippunfold_t2:
    """Export segmentation volume for a subject to TSV"""
    input:
        segs=expand(
            'hippunfold_highresT2/hippunfold/sub-{subject}/anat/sub-{subject}_hemi-{hemi}_space-cropT2w_desc-subfields_atlas-{atlas}_dseg.nii.gz',
                hemi=['L','R'], allow_missing=True),
        lookup_tsv='resources/hippunfold_desc-subfields_atlas-{atlas}_dseg.tsv'
    output:
        tsv='results/hippunfold-t2-{atlas}atlas/sub-{subject}_volumes.tsv'
    script:
        "scripts/gen_volume_tsv.py"




#get FS segmentation back into the space of the highresT2 (for comparing to others)
# need to apply the inverse of the T1_to_T2.lta 
rule convert_xfm_lta_itk:
    input:
        lta = '/home/ROBARTS/alik/graham/scratch/3_processing/hippunfold_HCA100/freesurfer_highresT2/hca100_freesurfer/freesurfer/sub-{subject}/mri/transforms/T1_to_T2.v21.lta',
#        lta = 'freesurfer/sub-{subject}/mri/transforms/T1_to_T2.v21.lta',
    output:
        itk = 'freesurfer/sub-{subject}/mri/transforms/T1_to_T2.v21.itk.txt'
    shell:
        'lta_convert --inlta {input.lta} --outitk {output.itk}'

rule transform_fs_to_t2:
    input:
        itk = 'freesurfer/sub-{subject}/mri/transforms/T1_to_T2.v21.itk.txt',
        seg = 'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.nii.gz',
        ref = 'hippunfold_highresT2/hippunfold/sub-{subject}/anat/sub-{subject}_desc-preproc_T2w.nii.gz'
    output:
        seg = 'freesurfer/sub-{subject}/mri/{hemi}.hippoAmygLabels-T2.v21.CA.space-T2w.nii.gz'
    shell:
        'antsApplyTransforms -i {input.seg} -r {input.ref} -t {input.itk} -o {output.seg} -n NearestNeighbor'
        

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



