function out =segMain(imOrig, configIn)
% Combined morphological-spectral multimodal segmentation
% Segments Image (imOrig) as a colour image if it is colour RGB image or in
% gray scale only if it is a gray image (one plane only).  
% 
% Inputs -  ImOrig: image to be segmented (should be in range 0-255)
%           configIn.Levels: no. of decomposition levels in the wavelet transform
%               (range: non-zero integers)
%           configIn.hmindepthfactor: controls the amount of the initial segmentation 
%               (range: 0->1)	
%           configIn.merge: perform merging(1) or not(0)
%           configIn.filterNoErode: do not know!
%
% Outputs - out.overlay: final segmentation overlayed on original image
%           out.intoverlay: initial segmentation (before merging) overlayed on 
%               original image
%           out.map: integer segmentation maps after merging (no borders)
%           out.intmap: initial integer segmentation maps before merging
%               (borders=0)
%           out.gsurf: gradient map

levels = configIn.levels;
hminfactor = configIn.hmindepthfactor;
merge = configIn.merge;
filterNoErode = configIn.filterNoErode;

imOrig = double(imOrig);

global globalmax;
globalmax = zeros(levels.*6,1); %store max value in each subband for each image

%% GET HIGHPASS SUBBANDS
[~,imFreqd] = dtwavexfm2(double(imOrig),levels,'antonini','qshift_06');

%% GET TEXTURE GRADIENT %%
% -------------------------------------------------------------------------
for n=1:levels
    imFreqd{n} = abs(imFreqd{n}(:,:,[1 6 3 4 2 5]));  %rearrange for backward compatibility with older code
    globalmax((n-1)*6+1:n*6) = squeeze(max(max(imFreqd{n},[],2),[],1))./(2.^n);
end

col = 0; %col=0 for all images; col=1 for CbCr channels; see 'texgradd.m'
[gsurf,imFreq] =segprotomm(imFreqd,double(imOrig), col, filterNoErode); %segment each image

clear imFreqd someIm;


%% SEGMENT %%
sed = strel('square',3);

% initial segmentation
gradmed = median(gsurf(:));
map_t = watershed(imhmin(gsurf,hminfactor*gradmed));
clear gradmed

%shave off nobbly bits in watersheds (trust me, I know what I'm talking
%about)
map2 = zeros(size(map_t));
for k=1:max(map_t(:))
    map2(imclose(map_t==k,sed)) = k;
end

[L,N] = superpixels(imOrig,500);
map2 = uint16(L); 
intmap = map2; %saved for display later



%% MERGE %%
if( merge )
    map = mergeRegions(map2,imFreq,imOrig,configIn); %merging regions
else
    map = map2;
end


%% POST PROCESS %%
overlay = uint8(imOrig);
intolay = uint8(imOrig);

overlayBand = overlay;
overlayBand(map==0) = 255;
overlay = overlayBand;
intolayBand = intolay;
intolayBand(intmap==0) = 255;
intolay = intolayBand;


%% This part is for removing region boudaries (map==0).

temp = map;
bw = ones(size(temp));
border = find(temp==0);
bw(border)=0;
[~, ind] = bwdist(bw);
temp(border) = temp(ind(border));
map = temp;

out.intolay = intolay;out.overlay = overlay;
out.intmap = intmap;out.map = map;
out.gsurf = gsurf;
