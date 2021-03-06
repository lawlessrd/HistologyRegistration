function [jpg_mask] = mask_hist(RGB_img,thresh,check)
% RGB img:
% threshold: theshold for keeping voxels in mask
%check: set check = 1 to view results at end of function


jpg_mask = zeros(size(RGB_img,[1,2,4]));

%generate circle in center of image
[m,n]=size(RGB_img(:,:,1,1));
[row,col]=meshgrid(1:n, 1:m);
radius = 200;

cir_mask = (row - n/2).^2/2^2 + (col - m/2).^2 ...
    <= radius.^2;

windowSize = 101;
kernel = ones(windowSize) / windowSize ^ 2;

for i = 1:size(RGB_img,4)
    lab_he = rgb2lab(RGB_img(:,:,:,i));
    ab = lab_he(:,:,2:3);
    ab = im2single(ab);
    nColors = 2;
    % repeat the clustering 3 times to avoid local minima
    pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
    mask1 = pixel_labels==1; mask2 = pixel_labels==2;
    
    %check which mask has more values in the center of image
    chck1 = sum(cir_mask.*mask1,'All');
    chck2 = sum(cir_mask.*mask2,'All');
    
    if chck1 > chck2
        S=1;
    else
        S=2;
    end

    % smooth binary mask to remove rough edges and fill holes
    tmp_mask=imfill(pixel_labels == S,'holes');
    
    blurryImage = conv2(single(tmp_mask), kernel, 'same');
    jpg_mask(:,:,i) = blurryImage > thresh; % Rethreshold
    
    
end

if check == 1
    for i = 1:size(RGB_img,4)
        figure(1); imshow(jpg_mask(:,:,i));
        title('Histology segmentation');
        pause(0.3);
    end
    %close(1);
end

