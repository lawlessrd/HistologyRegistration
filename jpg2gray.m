function [jpg_img,RGB_img] = jpg2gray()
%% Load in jpg images, export grayscale and RGB versions


% store jpg and convert to grayscale
jpg_files = dir('*.jpg');
jpg_size = size(imread(jpg_files(1).name));

jpg_img = zeros([jpg_size(1:2),length(jpg_files)]);
RGB_img = zeros([jpg_size,length(jpg_files)]);

for i = 1:length(jpg_files)
    [RGB, map, alpha] = imread(jpg_files(i).name);
    [rows, columns, numberOfColorChannels] = size(RGB);
    if size(RGB(:,:,1)) == size(jpg_img(:,:,i))
        jpg_img(:,:,i) = rgb2gray(RGB(:,:,:));
        RGB_img(:,:,:,i) = RGB(:,:,:);
    else
        jpg_img(:,:,i) = imresize(rgb2gray(RGB), size(jpg_img(:,:,i)));
        for j=1:3
            RGB_img(:,:,j,i) = imresize(RGB(:,:,j), size(RGB_img(:,:,j,i)));
        end
    end
end