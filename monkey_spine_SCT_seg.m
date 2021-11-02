function monkey_spine_SCT_seg(nii_file,subject)

%% load in nii

SCTfolder = sprintf('%s_SCT',subject);

if ~exist(SCTfolder, 'dir')
    mkdir(SCTfolder)
end

copyfile(nii_file,[SCTfolder, '/anat.nii']);
cd(SCTfolder);

unix('gzip anat.nii');

%% Begin SCT commands

% load in resampled nii
nii_info = niftiinfo('anat.nii.gz');
nii_img = imrotate(niftiread('anat.nii.gz'),-90);


% get centerline
% be sure to turn on manual mode in the viewer and select top and bottom slices
unix('sct_get_centerline -i anat.nii.gz -c t2 -method viewer');

% create mask for cropping
unix('sct_create_mask -i anat.nii.gz -p centerline,anat_centerline.nii.gz -size 12mm');

% remove area outside mask from image
unix('sct_maths -i anat.nii.gz -o anat_zero.nii.gz -mul mask_anat.nii.gz');
nii_img_zero = imrotate(niftiread('anat_zero.nii.gz'),-90);

% run segmentation
unix('sct_propseg -i anat_zero.nii.gz -c t1 -init-centerline anat_centerline.nii.gz -rescale 2'); 

nii_seg = imrotate(niftiread('anat_zero_seg.nii.gz'),-90);

