### MRI ###

# resample mri data
sct_resample -i SpineAnatImg.nii -vox 256x256x20

# get centerline (requires manual input)
# be sure to turn on manual mode in the viewer and select top and bottom slices
sct_get_centerline -i SpineAnatImg_r.nii -c t2 -method viewer

#fsleyes SpineAnatImg_r.nii -cm greyscale /Users/dylanlawless/Desktop/Histology/jpg2nii/SpineAnatImg_r_centerline.nii.gz -cm red -a 70.0 &

# create mask for cropping
sct_create_mask -i SpineAnatImg_r.nii -p centerline,SpineAnatImg_r_centerline.nii.gz -size 12mm

# remove area outside mask from image
sct_maths -i SpineAnatImg_r.nii -o SpineAnatImg_r_zero.nii -mul mask_SpineAnatImg_r.nii 

#fsleyes SpineAnatImg_r_zero.nii &

# propseg for spinal cord segmentation 
sct_propseg -i SpineAnatImg_r_zero.nii -c t1 -init-centerline SpineAnatImg_r_centerline.nii.gz -rescale 2 

# If the segmentation fails at some location (e.g. due to poor contrast between spinal cord and CSF), edit your anatomical image (e.g. with fslview) and manually enhance the contrast by adding bright values around the spinal cord for T2-weighted images (dark values for T1-weighted). Then, launch the segmentation again.
# https://sourceforge.net/p/spinalcordtoolbox/discussion/help/thread/346979d1/

### Histology ###
# be sure to turn on manual mode in the viewer and select top and bottom slices
sct_get_centerline -i hist_nii_r.nii -c t1 -method viewer

# GET SEG FROM MATLAB #
#sct_propseg -i hist_nii.nii -c t1 -init-centerline hist_nii_centerline.nii -rescale 2 -radius 2

sct_register_multimodal -i hist_nii_r.nii -iseg hist_nii_r_seg.nii \
  -d SpineAnatImg_r_zero.nii -dseg SpineAnatImg_r_zero_seg.nii \
  -o hist2mri.nii.gz -owarp warp_hist2mri.nii.gz \
  -param step=1,type=seg,algo=centermass,metric=MeanSquares:step=2,type=seg,algo=columnwise,metric=MeanSquares
