%upsample the texture subbands
% EG szmat = [size(tvolcell{1}(:,:,1)).*2 , 6.*levels];
% tvol = quickinterp(tvolcell,szmat);
% tvol = tvol(1:size(inimage,1),1:size(inimage,2),:);

function tvol = quickinterp(tvolcell,szmat,varargin);
%new version of the quick interpolation.  this can now handle arbitrary
%image dimensions, by choosing the padding and filter length (odd or even)
%to suit the desired final size
%RO'C 4/10/03


    tvol = zeros(szmat(1:2));
    for n=1:length(tvolcell);
        [M,N] = size(tvolcell{n}(:,:,1));
        zpadsize = ([M,N]-1).*(2.^n)+1;
        tvolpad = zeros([zpadsize(1),N,6]);
        
        %     inds1 = kron([1:M],ones(1,2^n));
        %     inds2 = kron([1:N],ones(1,2^n));
        %offset = 2^(n-1);
        start = (n-1)*6+1;
        tvolpad(1:2.^n:end,:,:)=tvolcell{n};
        filtsz = repmat(2^(n+2)+6,[1 2]);
        adjust = rem(szmat(1:2),2);
        filtsz = filtsz+adjust;
        Bc = fir1(filtsz(1)-1,1/(2^n)); %have made filters approx 2 times 
        Br = fir1(filtsz(2)-1,1/(2^n)); %longer than decimation factor
        Br = Br.*2^n; %fix gain because we have put in so many zeros
        Bc = Bc.*2^n;
        extnc = 0.5*(szmat(1)+filtsz(1)-zpadsize(1)-1);
        extnr = 0.5*(szmat(2)+filtsz(2)-zpadsize(2)-1);
        for m=1:6;
            %extend symmetrically about the last pixel
            temp = zeros(szmat(1),zpadsize(2));
            temp(:,1:2.^n:zpadsize(2)) = conv2(tvolpad([ [extnc+1:-1:2],[1:end],[end-1:-1:end-extnc] ],:,m),Bc','valid');
            tmptvol = conv2(temp(:,[ [extnr+1:-1:2],[1:end],[end-1:-1:end-extnr] ]),Br,'valid');
            tmptvol = max(tmptvol,0);
            tvol = tvol + tmptvol;
        end;
    end;
    
