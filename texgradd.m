function outim =texgradd(tvol, col)
szmat = size(tvol{1}(:,:,1))*2;
szmat = [szmat, length(tvol).*6];

gvol = tvol;

for k=1:length(tvol);
    for n=1:6;
        if( col == 0 )
            temp= gradim(tvol{k}(:,:,n),1.5); % Y channel do 1.5, 3 for others
        else
            temp= gradim(tvol{k}(:,:,n),3); % Y channel do 1.5, 3 for others
        end
        norm = max(temp(:));
        temp = temp./norm;
        norm = sum(temp(:).^2)./prod(size(tvol{k}(:,:,n)));
        temp = temp./norm;
        
        gvol{k}(:,:,n) = temp;%./2.^k;
    end;
end;

outim = quickinterpPRH(gvol,szmat);
% Now, clean up, canny style
% Using pixels over a  high threshold as the marker image, use morph.
% reconstruction to include lower, but connected gradient peaks
omed = median(outim(:));
if omed>0;
    outim = outim./omed;
end;
marker = outim;
%        marker(marker<(0.1)) = 0;
marker(marker<1) = 0;        
outim=imreconstruct(marker,outim,4);
% Now floor values beneath the low threshold
%        temp(temp<0.05)=0.05;
outim(outim<0.2)=0;
