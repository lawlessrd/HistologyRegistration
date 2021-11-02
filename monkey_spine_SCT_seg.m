function [nii_info,nii_img,nii_img_zero,nii_seg] = monkey_spine_SCT_seg(subject)

%% load in nii
% if ~exist('anat_r.nii.gz','file')
%     cd(nii_dir);
% end

SCTfolder = sprintf('%s_SCT',subject);

if ~exist(SCTfolder, 'dir')
    mkdir(SCTfolder)
end


%% Begin SCT commands

cd(SCTfolder);

% load in resampled nii
nii_info = niftiinfo('anat_r.nii.gz');
nii_img = imrotate(niftiread('anat_r.nii.gz'),-90);


% get centerline
% be sure to turn on manual mode in the viewer and select top and bottom slices
unix('sct_get_centerline -i anat_r.nii.gz -c t2 -method viewer');

% create mask for cropping
unix('sct_create_mask -i anat_r.nii.gz -p centerline,anat_r_centerline.nii.gz -size 12mm');

% remove area outside mask from image
unix('sct_maths -i anat_r.nii.gz -o anat_r_zero.nii.gz -mul mask_anat_r.nii.gz');
nii_img_zero = imrotate(niftiread('anat_r_zero.nii.gz'),-90);

% run segmentation
unix('sct_propseg -i anat_r_zero.nii.gz -c t1 -init-centerline anat_r_centerline.nii.gz -rescale 2'); 

nii_seg = imrotate(niftiread('anat_r_zero_seg.nii.gz'),-90);

