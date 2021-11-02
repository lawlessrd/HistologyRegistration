%% To do: 
% test masking on different histology types (make into function)
% find appropriate C2 staining 
% get code to generate niftis from raw data
% anisotropic voxels
% fill holes/dilate hist mask

clear all; 
close all;
%% Set paths and directories
Path = getenv('PATH');
if ~contains(Path,'/Users/dylanlawless/sct_5.3.0/bin')
    setenv('PATH', [Path ':/Users/dylanlawless/sct_5.3.0/bin']);  %Edit rdl 4/19/2018 getenv was just PATH
end

home = pwd;
nii_dir = 'sqm6356';
jpg_dir = '6356 C2';
subject = '6356';

nii_file = 'SpineAnatImg.nii';

%% Use SCT to resample and segment MRI

[nii_info,nii_img,nii_zero,nii_seg] = monkey_spine_SCT_seg(nii_dir,nii_file,subject);


%% Load in jpg images, export grayscale and RGB versions
% remove unwanted jpgs from directory before running
cd(home);

[jpg_img,RGB_img] = jpg2gray(jpg_dir);
   
%% manipulate jpg to better resemble mri

% inv colors to better represent MRI data
jpg_img_inv = imcomplement(jpg_img);
jpg_img_s = zeros([size(jpg_img,[1,2]),size(nii_zero,3)]);

% set range of values to match nii
% need to edit for image spacing here, will only use jpgs up to max number
% of nifti slices
for i = 1:size(nii_zero,3)
    nii_tmp = nii_zero(:,:,i);
    jpg_img_s(:,:,i) = double(max(nii_tmp(:)) * mat2gray(jpg_img_inv(:,:,i)));
end

%% mask jpg 

thresh = 0.5;
check = 1;

[jpg_mask] = mask_hist(RGB_img,thresh,check);


%% shrink hist and zero pad to match size of MRI

fake_nii = zeros(size(nii_img));
mask_nii = zeros(size(nii_img));

for i = 1:size(nii_img,3)
    
    [row,col]=find(nii_seg(:,:,i) ==1);
    
    row_min = min(row);
    row_max = max(row);
    row_diff = row_max-row_min;
    
    col_min = min(col);
    col_max = max(col);
    col_diff = col_max-col_min;
    
    jpg_shrink{i} = imresize(jpg_img_s(:,:,i), [row_diff+1,col_diff+1]);
    jpg_mask_shrink{i}=imresize(jpg_mask(:,:,i), [row_diff+1,col_diff+1]);
    
    jpg_mask_shrink{i}(jpg_mask_shrink{i} > 0.3) = 1;
    jpg_mask_shrink{i}(jpg_mask_shrink{i} < 0.3) = 0;
    
    %zero pad
    jpg_mask_zero = zeros(size(nii_img(:,:,1)));
    jpg_mask_zero(row_min:row_max,col_min:col_max) = jpg_mask_shrink{i};
    
    jpg_shrink_zero = zeros(size(nii_img(:,:,1)));
    jpg_shrink_zero(row_min:row_max,col_min:col_max) = jpg_shrink{i};
    
    nii_tmp = zeros(size(nii_img(:,:,i))); %nii_tmp = nii_img(:,:,i);
    mask_tmp = zeros(size(nii_img(:,:,i)));
    
    fake_nii(:,:,i) = jpg_shrink_zero;
    mask_nii(:,:,i)=jpg_mask_zero;
end

%% Check results and save new nii

figure(); subplot(121); dispimg(nii_img(:,:,1));
subplot(122); dispimg(fake_nii(:,:,1));

cd(home);
cd(nii_dir);
cd(dir('*_SCT').name)

niftiwrite(single(imrotate(fake_nii,90)), 'hist_nii_r.nii',nii_info);
niftiwrite(single(imrotate(mask_nii,90)), 'hist_nii_r_seg.nii',nii_info);

%% Register images using SCT
inp = sprintf(['sct_register_multimodal -i hist_nii_r.nii -iseg hist_nii_r_seg.nii ' ...
  '-d anat_r_zero.nii -dseg anat_r_zero_seg.nii ' ...
  '-o hist2mri.nii.gz -owarp warp_hist2mri.nii.gz ' ...
  '-param step=1,type=seg,algo=centermass,metric=MeanSquares:' ...
  'step=2,type=seg,algo=columnwise,metric=MeanSquares']);

unix(inp);

%% Check results
% hist_nii size is a bit off
% correct spacing

%% adjust contrast of cropped nii to improve segmentation
% nii_crop = niftiread('SpineAnatImg.nii');
% nii_crop_info = niftiinfo('SpineAnatImg.nii');
% 
% nii_cont = zeros(size(nii_crop));
% 
% for i = 1:size(nii_crop,3)
% nii_tmp=nii_crop(:,:,i);
% 
% % contrast
% meanV = mean2(nii_tmp);
% newV = meanV * 2 + 5 * (nii_tmp - meanV); % Increase contrast by factor of 1.5, birghtness by factor of 2
% 
% nii_cont(:,:,i) = newV;
% end
% 
% nii_crop_info.Filename = '/Users/dylanlawless/Desktop/Histology/jpg2nii/sqm6356/SpineAnatImg_cont.nii';
% niftiwrite(single(nii_cont),'SpineAnatImg_cont.nii',nii_crop_info);
% 
% figure(); subplot(121); imshow(nii_tmp,[min(nii_tmp(:)),max(nii_tmp(:))]); colormap gray;
% subplot(122); imshow(newV,[min(nii_tmp(:)),max(nii_tmp(:))]);

%% manually mask spinal cord image

% nii = niftiread('SpineAnatImg_cont_r.nii');
% nii_info = niftiinfo('SpineAnatImg_cont_r.nii');
% 
% for i = 1:size(nii,3)
%     figure(1);
%     pause(0.00001);
%     frame_h = get(handle(gcf),'JavaFrame');
%     set(frame_h,'Maximized',1);
%     dispimg(nii(:,:,i));
%     title('SC')
%     sc_mask(:,:,i) = roipoly;
%     %title('CSF');
%     %csf_mask(:,:,i) = roipoly;
% end
% 
% nii_info.Filename = '/Users/dylanlawless/Desktop/Histology/jpg2nii/sqm6356/nii_seg.nii';
% niftiwrite(single(sc_mask), 'nii_seg.nii',nii_info);