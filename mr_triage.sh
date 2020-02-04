#!/bin/sh
Usage() {
    cat <<EOF

MRtriage based on symmetry comparison between both hemispheres v0.1 (alpha)
To run run it as as a daemon install fswatch (brew install fswatch)
Developed by gABoCas 16.02.2016


Usage: mr_triage <input_folder/file> <bet_flag>

e.g.   mr_triage ~/dcm_dir 1

EOF
    exit 1
}

BASEDIR="$1"
[ "$1" = "" ] && Usage
# Global vars def
dcm2nii_path=/Applications/mricro
mni_brain=${FSLDIR}/data/standard/MNI152_T1_2mm_brain

mcverter -f fsl -d -n  -F +SeriesDescription -o ${BASEDIR} ${BASEDIR}
nii_file=`ls ${BASEDIR}/*.nii`

in_dir=`dirname $nii_file`
nii_file=`remove_ext $nii_file`
# If the bet flag is set to 1, the nifti file is brain extracted (the input for flirt should be brain extracted)
bet_flag=$2; shift
echo ${bet_flag}
if [ $bet_flag = 1 ] ; then
	bet $nii_file ${nii_file}_brain -R -f 0.5 -g 0 -m
	nii_file=${nii_file}_brain
	mask_file=${nii_file}_mask
	fslmaths $mask_file -eroF $mask_file
	csf_thr=`fslstats $nii_file -k $mask_file -s`
	fslmaths $nii_file -mas $mask_file -uthr $csf_thr -dilM -bin ${nii_file}_csf_mask
	fslmaths $mask_file -mas $mask_file -sub ${nii_file}_csf_mask -bin $mask_file
	wait
fi
#spatial smoothing using a Gauss filter FWHM=2.355 sigma
fslmaths $nii_file -s 1.3 $nii_file 
#Registration to standard space
flirt -in $nii_file -ref $mni_brain -out ${nii_file}_mni-2mm -omat ${nii_file}_mni-2mm.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear
# Apply transformation to the mask
flirt -in $mask_file -ref $mni_brain -out $mask_file -init ${nii_file}_mni-2mm.mat -applyxfm
nii_file=${nii_file}_mni-2mm

# Splitting of both hemispheres 
fslroi ${nii_file} ${nii_file}_r 0 45 0 -1 0 -1
fslroi ${nii_file} ${nii_file}_l 45 45 0 -1 0 -1
fslroi $mask_file ${mask_file}_r 0 45 0 -1 0 -1
fslroi $mask_file ${mask_file}_l 45 45 0 -1 0 -1
fslswapdim ${nii_file}_l -x y z ${nii_file}_l_flipped
fslswapdim ${mask_file}_l -x y z ${mask_file}_l_flipped
fslmaths ${mask_file}_r -mul ${mask_file}_l_flipped -bin ${mask_file}_no_csf

# NCC
./ncc.sh NCC_r-fl ${nii_file}_r ${nii_file}_l_flipped ${mask_file}_no_csf
./ncc.sh NCC_r ${in_dir}/NCC_r-fl ${nii_file}_r ${mask_file}_no_csf
ncc_r=`fslstats ${in_dir}/NCC_r -R | awk '{print $2}'`
./ncc.sh NCC_fl ${in_dir}/NCC_r-fl ${nii_file}_l_flipped ${mask_file}_no_csf
ncc_fl=`fslstats ${in_dir}/NCC_fl -R | awk '{print $2}'`
ncc_thr=`echo "(${ncc_r}+${ncc_fl})/6" | bc -l`
fslswapdim ${in_dir}/NCC_fl -x y z ${in_dir}/NCC_l
fslroi ${nii_file} ${nii_file}_line 90 1 0 -1 0 -1
fslmerge -x ${in_dir}/NCC_map ${in_dir}/NCC_r ${in_dir}/NCC_l ${nii_file}_line
#new
ncc_map_nt=${in_dir}/NCC_map_nt
cp ${in_dir}/NCC_map.nii.gz ${ncc_map_nt}.nii.gz
fslmaths ${in_dir}/NCC_map -thr $ncc_thr -eroF -fmedian -dilM -bin ${in_dir}/NCC_map
##
zero_ncc=`fslstats ${in_dir}/NCC_map -R | awk '{print $2}' | bc -l`
if [ $zero_ncc == 0 ];then
	echo "performing alternative post-processing"
	fslmerge -x ${in_dir}/NCC_map ${in_dir}/NCC_r ${in_dir}/NCC_l ${nii_file}_line
	fslmaths ${in_dir}/NCC_map -thr $ncc_thr -fmedian -dilD -eroF -fmedian -bin ${in_dir}/NCC_map
fi
#new
fslmaths ${ncc_map_nt} -uthrP 1 -eroF -bin ${ncc_map_nt}_u1
fslmaths ${ncc_map_nt} -thrP 75 -bin ${ncc_map_nt}_75
fslmaths ${mask_file} -thr 0.05 -bin -add -1 ${mask_file}_neg_bk
fslmaths ${ncc_map_nt}_u1 -add 1 -thr 1.1 -add ${in_dir}/NCC_map -mas ${mask_file} -add ${mask_file}_neg_bk  ${in_dir}/rw_markers
##
rm ${nii_file}_line*
#fslmaths ${in_dir}/NCC_map -thr $ncc_thr -eroF -fmedian -dilM -bin ${in_dir}/NCC_map
#echo $ncc
