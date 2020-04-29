function testseg
%   Author: Paul Hill, Rob O'Callaghan, Pui
%           VILab
%           the University of Bristol, UK
%
%   Copyright (c) 2002--2020 by  Paul Hill, Rob O'Callaghan, Pui
%
%   References:
%    [1] R. O'Callaghan,  and D.R. Bull. "Combined morphological-spectral 
%        unsupervised image segmentation." IEEE transactions on image
%        processing 14, no. 1 (2004): 49-62.
%    [2] Hill, P. R., Canagarajah, C. N., & Bull, D. R. (2003). Image 
%        segmentation using a texture gradient based watershed transform.
%        IEEE Transactions on Image Processing, 12(12), 1618-1633.
%
% Function that is the top level test for UoB segmentation code
%
% USAGE:
%   testSeg
% INPUT:
%   -
% OUTPUT:
%   out.gurf is the gradient surface
%   out.segMap is the segmentation map (one number per region)
%   out.intmap is the pre merged segemntation map
%   out.overlay is the segementation edges overlayed on the Image
%   out.intolavy is the same but onthe pre-merged segmentation map


%% INPUT IMAGE MUST HAVE DIMENSIONS OF A POWER OF 2 CURRENT (i.e. 128x128, 256x256)
%% ALSO number of regions are limited to under ~1000 due to memory restrictions

addpath('./dtcwt/');
Im_RGB = imread('lena.png');
if length(size(Im_RGB)) == 3
	Im = rgb2ycbcr(Im_RGB);
	Im = Im (:,:,1);
else
	Im = Im_RGB;
end

configIn.hmindepthfactor = 0.15;  
configIn.levels = 4;
configIn.t1= 0.9; configIn.t2= 0.8; 
configIn.merge = 1;
configIn.nonErode = 0;
configIn.filterNoErode = 0;

%t1 and t2 are usual merging thresholds
%levels are the number of wavelet levels
%merge is whether to do the merging of the initial segmentation
%hmindepthfactor controls the amount of the initial segmentation it can be set between 0 and 1 

out = segMain(Im, configIn);

out.overlay; out.intolay;
out.map;  out.intmap;
out.gsurf;

imagesc(uint8(out.map));



