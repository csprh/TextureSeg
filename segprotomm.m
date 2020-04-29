function[gradsurf,tvol]=segprotomm(tvolcell,inimage,col,filterNoErode)


% if directional order-statistic fitering isn't applied to non erosion
% version, store the original tvolcell before filtering process.
if ~filterNoErode
    % preparing for non erosion version
    tvolcellNoErodeCell = tvolcell;
end

% applying directional order-statistic filtering
tvolcell = texmedfiltd(tvolcell);
if filterNoErode
    tvolcellNoErodeCell = tvolcell;
end

levels = length(tvolcell);
tgradim = texgradd(tvolcell, col);
tgradim = tgradim(1:size(inimage,1),1:size(inimage,2));

% erode back texture response so as not to stifle valid gradients at edges
% of texture regions

se = strel('square',3);

for n=1:levels;
    
    % <PUI> Not applying erosion so that details aren't washed out
    tvolcellNoErodeCell{n} = tvolcellNoErodeCell{n}/2^n;
    
    for m=1:6;
        tvolcell{n}(:,:,m)=imerode(tvolcell{n}(:,:,m),se)./2^n;
    end;
end;

szmat = [size(tvolcell{1}(:,:,1)).*2,6.*levels];

% <PUI> will use this in merging region
if ~filterNoErode
    tvol = quickinterp(tvolcell,szmat);
else
    tvol = quickinterp(tvolcellNoErodeCell,szmat);
end
tvol = tvol(1:size(inimage,1),1:size(inimage,2),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tvolcell;

texen = gettexen(tvol);

if col == 0
    greygradim = gradim(inimage,1.5);
else
    greygradim = gradim(inimage,3);
end;

gmed = median(greygradim(:));
greygradim(greygradim<0.5*gmed)=0;  % to kill off noise before we modify and potentially amplify it

if gmed>0;
    greygradim = greygradim./gmed;
end;

tx2 = (texen/2-7); %1.5-3, 2-7

clear texen;
tx2 = exp(max(tx2,0));
greygradim = greygradim./tx2;  
clear tx2;

% remove anything below a threshold - this is not noise removal, but low
% amplitude basin removal for watershed, hence we clamp to the threshold not
% zero
tgradim(tgradim<(0.1*max(tgradim(:)))) = 0.1*max(tgradim(:)); 

gradsurf = tgradim.*4+greygradim;%./absolutetexfac;     mutliplied by four
%gradsurf = tgradim/max(tgradim(:)) + greygradim/max(greygradim(:)); %added by Essa to make equal contribution of both texture and gradient
%gradsurf = tgradim; % <Essa> only gradient and no texture


clear tgradim greygradim;


% this is orignal function
gradsurf = filter2(ones(3),gradsurf);

