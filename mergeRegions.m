function outMap = mergeRegions(map,inFreq,inIm,configIn)
% Combined morphological-spectral multimodal segmentation
% Segments Image (imOrig) as a colour image if it is colour RGB image or in
% gray scale only if it is a gray image (one plane only).  
% 
% Inputs -  
%           map: segmentation map
%           inFreq: frequency content (w,h,numberOfChannels)
%           inIm: image to be segmented (should be in range 0-255)
%           configIn.t1 and configIn.t2: are usual merging thresholds (range: 0->1)
% Outputs - outMap: merged segmentation map

t1 = configIn.t1;
t2 = configIn.t2;

binsperdim=64;

origmap = map;
regions=nonzeros(unique(map(:)));
numreg = length(regions);

adjacemat = spalloc(numreg,numreg,3*numreg);   %Sparse version
dists1 = adjacemat;

se1 = strel('line',6,0);
se2 = strel('line',6,90);
se3 = strel('square',3);
se4 = strel('square',5);

%% Erode Regions so we get a central region to characterise
%generate eroded region map, which is only used to train the classes for
%region difference calculation
erodemap = zeros(size(map)+2);
map4erode = erodemap;
map4erode(2:end-1,2:end-1) = map;  %need to pad coz otherwise regions on periphery will erode towards edge!

for k=1:numreg
    template=(map4erode==k);
    prevtemplate=template;
    numerodes=0;
    while((nnz(template)>(6*binsperdim))&(numerodes<7)) 
        prevtemplate=template;
        template=imerode(template,se4);
        numerodes=numerodes+1;
    end
    erodemap(prevtemplate)=k;
end
clear map4erode;
erodemap = erodemap(2:end-1,2:end-1);


keepinds = find(erodemap);

featFreq = reshape(inFreq,size(inFreq,1)*size(inFreq,2),size(inFreq,3));
featIm = reshape(inIm,size(inIm,1)*size(inIm,2),size(inIm,3));
clear inFreq inIm;


%% Quantize texture feature set  (involves scaling into uniform range)
%make the bins equal width, over the range of the image

maxvec = max(featFreq,[],1);
global globalmax;
maxvec = max(globalmax(:).'./2,maxvec);

binwidths = maxvec./binsperdim;
%featFreq = floor(featFreq./repmat(binwidths,size(featFreq,1),1))+1;
for k=1:size(featFreq,1)
    featFreq(k,:) = floor(featFreq(k,:)./binwidths)+1;
end
featFreq(featFreq>binsperdim)=binsperdim;  %just to force maxima into the bins


tempFeats = zeros(size(featIm));
featIm = featIm(keepinds,:);

maxvec = 255.*ones(1,size(featIm,2));

binwidths = maxvec./binsperdim;
featIm = floor(featIm./repmat(binwidths,size(featIm,1),1))+1;
featIm(featIm>binsperdim)=binsperdim;  %just to force maxima into the bins

tempFeats(keepinds,:) = featIm;
featIm = tempFeats;
clear tempFeats;


%% Find Adjacency Matrix (Dilate in 2D and Find Overlap)
for k=1:numreg
    overlap = imdilate((map==regions(k)),se1); %for watershed pixels
    overlap = overlap|imdilate((map==regions(k)),se2);
    map(map==k)=0;
    adjvec = map(overlap);
    adjvecun=unique(adjvec(adjvec>0));
    for i=1:size(adjvecun)
        adjacemat(adjvecun(i),k)=nnz(adjvec==adjvecun(i));
    end
end
map=origmap;
clear origmap;
adjacemat = adjacemat + adjacemat';

%% Find calculate "distances" between adjacent regions
[xpos,ypos] = find(triu(adjacemat));
for k=1:length(xpos)

    template1=(erodemap==xpos(k));
    template2=(erodemap==ypos(k));
    reg1t = featFreq(template1(:),:);
    reg2t = featFreq(template2(:),:);
    reg1g = featIm(template1(:),:);
    reg2g = featIm(template2(:),:);
    
    reg1t(reg1t==0) = 1;    reg2t(reg2t==0) = 1;
    reg1g(reg1g==0) = 1;    reg2g(reg2g==0) = 1;
    if ((sum(sum(reg1t==0))+sum(sum(reg2t==0))+sum(sum(reg1g==0))+sum(sum(reg2g==0)))>0)
        dists1(xpos(k),ypos(k)) = 1000000000000;
    else
        
        % <PUI> using matlab version - added weight to texdiff according to
        % mean magnitude of each subband
        [texdiff,greydiff]=distcalcmhMAT(reg1t,reg2t,reg1g,reg2g,binsperdim);
             
        % <Essa> measure the regularity of boundary using fractal dimension
        regdiff = calRegIndex(map, xpos, ypos, k);
        
        dists1(xpos(k),ypos(k)) = max([greydiff,texdiff,regdiff]);
    end

end
clear featFreq featIm erodemap;

dists = dists1+dists1.';

aff = dists;


[ai, aj] = find(aff);
for yaff = 1: size(ai,1)
        aff(ai(yaff),aj(yaff)) = 1 - aff(ai(yaff),aj(yaff));
end

saff = size(aff);
for yaff = 1: saff(1)
        aff(yaff,yaff) = 1;
end

inInds = 1:saff(1);
global globalInds;
globalInds = 1:saff(1);

% <PUI> try without concerning how much two areas connect - seem to produce
% more number of regions (fewer merged regions)
% specsplit5(aff,adjacemat,t1,t2, inInds);
specsplit5(aff.*adjacemat,adjacemat,t1,t2, inInds);

W = aff;


outMap = zeros(size(map));
uniqueInds = unique(globalInds);
for k=1:length(uniqueInds)
    fInds = find(globalInds==uniqueInds(k));
    for n=1:length(fInds)
        outMap(map==fInds(n))= k;
    end
end

for k=0:max(outMap(:))
    outMap(imclose(outMap==k,se3))=k;
end
    
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [texturediff,intensediff] = distcalcmhMAT(r1t,r2t,r1g,r2g,binsperdim) %marginal histograms with either EMD or L1 measures


%texture first
r1 = r1t; r2 = r2t;
N=size(r1,1);
M=size(r2,1);
numdims = size(r1,2);

histObj1 = zeros(binsperdim,numdims);
histObj2 = histObj1;

for n=1:N
    for k=1:numdims
        histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
    end
end
for n=1:M
    for k=1:numdims
        histObj2(r2(n,k),k) = histObj2(r2(n,k),k) +1;
    end
end
histObj1 = histObj1./N;
histObj2 = histObj2./M;
%1-D EMD dist now

histObj1 = cumsum(histObj1,1);
histObj2 = cumsum(histObj2,1);

histObj1 = abs(histObj1-histObj2);
dvec = sum(histObj1,1)./(binsperdim-1);

% <PUI> find weight from binary distance ----------------------------------
meanMag1 = mean(r1,1);
meanMag2 = mean(r2,1);
% binary pattern
binaryPattern1 = meanMag1(2:end) > meanMag1(1:end-1);
binaryPattern2 = meanMag2(2:end) > meanMag2(1:end-1);
% weight
weight = sum(binaryPattern1~=binaryPattern2);
weight = weight/length(meanMag1)*6;
% -------------------------------------------------------------------------

% final texture distance value
texturediff = min(max(dvec).*weight,1);

%now greyscale
r1 = r1g; r2 = r2g;
N=size(r1,1);
M=size(r2,1);
numdims = size(r1,2);

histObj1 = zeros(binsperdim,numdims);
histObj2 = histObj1;

for n=1:N
    for k=1:numdims
            histObj1(r1(n,k),k) = histObj1(r1(n,k),k) +1;
    end
end
for n=1:M
    for k=1:numdims
            histObj2(r2(n,k),k) = histObj2(r2(n,k),k) +1;
    end
end

histObj1 = histObj1./N;
histObj2 = histObj2./M;
%1-D EMD dist now

histObj1 = cumsum(histObj1,1);
histObj2 = cumsum(histObj2,1);

histObj1 = abs(histObj1-histObj2);
dvec = sum(histObj1,1)./(binsperdim-1);
intensediff = min(max(dvec).*2,1);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function regdiff = calRegIndex(map, xpos, ypos, k)

se1 = strel('square',3);
regionX = imdilate((map==xpos(k)),se1); %dilate first region
regionY = imdilate((map==ypos(k)),se1); %dilate second region
border = (map==0); %extract watershed lines
%the shared border is the intersection of previous regions
sharedBorder = regionX & regionY & border; 


%cropping the border region from the image
[r, c] = find(sharedBorder);
sharedBorder = sharedBorder(min(r):max(r), min(c):max(c));

%only if the boundry is long enough calculate the fractal dimension
if (max(r)-min(r))>5 || (max(r)-min(r))>5
    addpath('boxcount'); %include fractal library
    [n, r] = boxcount(sharedBorder);
    df = -diff(log(n))./diff(log(r));
    %disp(['Fractal dimension, Df = ' num2str(mean(df(:))) ' +/- ' num2str(std(df(:)))]);

    if mean(df(:))==1
        regdiff = 1;
    else
        regdiff = 0;
    end
    
else
    regdiff = 0;
end

return;