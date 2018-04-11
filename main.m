%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Implementation of 'Multispectral Demosaicking with Novel Guide Image 
%   Generation and Residual Interpolation'
%
%   Reference
%     Yusuke Monno, Daisuke Kiku, Sunao Kikuchi, Masayuki Tanaka,
%     and Masatoshi Okutomi,
%     'Multispectral Demosaicking with Novel Guide Image Generation and 
%     Residual Interpolation,'
%     Proc. of IEEE International Conference on Image Processing(ICIP2014),
%     pp.645-649, October, 2014.
%  
%   License 
%     This code is available only for reserch purpose.
%     You may not re-publish, re-transmit, or otherwise distribute directly 
%     or via links any data to any third person except for your personal, 
%     non-commercial use.
% 
%   Used array pattern  
%      G  R  G  Or G
%      B  G  Cy G  B
%      G  Or G  R  G
%      Cy G  B  G  Cy
%      G  R  G  Or G
%
% -------------------------------------------------------------------------
%   Copyright (C) 2014 Yusuke Monno. All rights reserved.
%   ymonno@ok.ctrl.titech.ac.jp
%
%   Dec. 18, 2014.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **** modified a little bit to implement: Single-Sensor RGB and NIR Image 
% Acquisition: Toward Optimal Performance by Taking Account of CFA Pattern, Demosaicking, and Color Correction
% by xavysp
% ****
clear all;
is_image = true; % true for read image and false for dataset (h5 file)

if is_image
    srcN = 1280;
    srcM = 720;
    path = 'TE3-RGBNIR-16_06-608.raw';
    f = fopen(path, 'r');
    I = fread(f, [srcN, srcM], 'uint16');
    fclose(f);
    I = I';
    I = double(I/(2^10-1));
    B = I(1:2:end, 1:2:end);
    IR = I(1:2:end, 2:2:end);
    G = I(2:2:end, 1:2:end);
    R = I(2:2:end, 2:2:end);
    vn(:,:,1)=R;
    vn(:,:,2)=G;
    vn(:,:,3)=B;
    vn(:,:,4)=IR;
    clear G R B IR;
    
    
else
    path = 'OMSIV_test_192.h5';
    data = h5read(path, '/data');
    label = h5read(path,'/label');
    label = permute(label,[4,3,2,1]);
    data = permute(data,[4,3,2,1]);
    disp(size(data));
    disp(size(label));
    vn = data(54,:,:,:);
    vn = reshape(vn,[192,192,4]);
    v = label(54,:,:,:);
    v = reshape(v,[192,192,3]);
    clear data label
end
% Load input images

% Degamma correction
gam = 2.2;
% Mosaicking and mask
% G R G O ..
% B G C G ..
% G O G R ..
% C G B G ..
% : : ::
pattern = 'bgnr';
pattern1 = 'bgbgngbggngrgrgrngbgbgbggrgrgrgnbgbgbgnggrgngrgrbgngbgbggrgrgrgr';
pattern2 = 'rgbggngnbgrggngn';
which_one = 1; % Choice 0 for pattern, 1 for pattern1 and 2 for pattern2

[mosaic mask] = mosaic_4band(vn,pattern1,which_one);
% subsampled each pixels  pag. 6 fig 3c

% Demosaicking
[green, rgb,nir] = MSIdemosaicking(mosaic,mask, which_one);

% Transform to sRGB image
% demsRGB = ROGCB2sRGB(demRGBOC(:,:,1),demRGBOC(:,:,4),demRGBOC(:,:,2),demRGBOC(:,:,5),demRGBOC(:,:,3));
% rgb = RGB2sRGB(rgb(:,:,1),rgb(:,:,2),rgb(:,:,3));

% Ggamma correction
rgb = norm_image(rgb,1);
nir = norm_image(nir,1);
rgb = gamma_correction(rgb,gam);
nir = gamma_correction(nir, gam);
rgb = norm_image(rgb,255);
nir = norm_image(nir,255);
rgb = uint8(rgb);
nir= uint8(nir);
rgb1 = imresize(rgb,size(I)/2);
% imshow(rgb);
% figure;
% imshow(nir)
imwrite(rgb,'visible.png');
imwrite(nir,'nir.png');

% Save result images
% imwrite(demRGBOC(:,:,1),sprintf('%s_demR.bmp',imagename));
% imwrite(demsRGB,sprintf('%s_demsRGB.bmp',imagename));

