function match_fov(home,nii_file,subject, varargin)
%% Check varargin for slice info
% varargin 1: starting nifti slice
% varargin 2: number of histology slices
% varargin 3: histology slice thickness (mm) (thickness+spacing)
if size(varargin,2) == 0
    fprintf('No slice info given, setting default values.');
    nii_slice = 1;
    hist_n_slices = 25;
    hist_thick = 0.3;
elseif size(varargin,2) == 1
    nii_slice = varargin{1};
    hist_n_slices = 25;
    hist_thick = 0.3;
elseif size(varargin,2) == 2
    nii_slice = varargin{1};
    hist_n_slices = varargin{2};
    hist_thick = 0.3;
elseif size(varargin,2) == 3
    nii_slice = varargin{1};
    hist_n_slices = varargin{2};
    hist_thick = varargin{3};
else
    fprintf('Exceeded number of inputs. Please repeat.');
end

%% load in nii

SCTfolder = sprintf('%s_SCT',subject);

if ~exist(SCTfolder, 'dir')
    mkdir(SCTfolder)
end

copyfile(nii_file,SCTfolder);
cd(SCTfolder);

%% Split and merge nii using fslsplit/fslmerge
nii_info = niftiinfo(nii_file);

nii_thick = nii_info.ImageSize(3)/nii_info.PixelDimensions(3);

num_slices = round(hist_n_slices/(nii_thick/hist_thick));

unix(sprintf('fslsplit %s -z',nii_file));

split_nii = dir('vol*.nii.gz');

merge_files = '';
for i = nii_slice:nii_slice+num_slices
    merge_files = append(merge_files, split_nii(i).name,' ');
end

unix(sprintf('fslmerge -z anat_slices%ito%i.nii.gz %s',nii_slice,nii_slice+num_slices,merge_files));

delete(split_nii.name);

% resample nii with sct_resample
unix(sprintf('sct_resample -i anat_slices%ito%i.nii.gz -vox %ix%ix%i -o anat_r.nii.gz', ...
    nii_slice,nii_slice+num_slices,nii_info.ImageSize(1) * 2,nii_info.ImageSize(2) * 2,hist_n_slices));

cd('..');

