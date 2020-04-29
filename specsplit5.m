function specsplit5(Win, Ain, t1, t2, inIndsin)
%now using new formulation for maximised mean association
%and with propagation of external edges for decision making
global globalInds;

queueitem{1} = Win; queueitem{2} = Ain;queueitem{3} = inIndsin;
queue{1} = queueitem;


% take item off back of queue and process in loop
%% This is much implemented as a "non recursive function" i.e.  as a Queue
while size(queue,2) > 0
    sizeOfQueue = size(queue,2);
    queueitem1 = queue{sizeOfQueue};
    W = queueitem1{1};
    A = queueitem1{2};
    inInds = queueitem1{3};
    queue = queue(1:sizeOfQueue-1);
    
    numreg=size(W,1);
    thresh1 = t1;
    thresh2 = t2;
    
    sw = size(W);
    for yw = 1: sw(1)
        W(yw,yw) = 0;
    end
    
    
    if (numreg>2);  %non-trivial cut computation
        
        d=sum(W,2);
        
        D = diag(d);
        S = diag(sum(A,2));
        
        Nr = (D-W);
        Dr = (S-A);
        [n,n] = size(Nr);
        if (n < 500)
            [V,evals]=eig(full(Nr),full(Dr)); %seems faster than using eigs to get just 2 eigenvectors
            [minval,minind] = min(abs(diag(evals)));
            if( (max(V(:,minind))-min(V(:,minind)))<1e-3 );
                evals(minind,minind) = inf;
                [minval,minind] = min(abs(diag(evals)));
            end;
            [Y,inds] = sort(V(:,minind));
        else
            [V,evals] = lobpcg(rand(n,5),Nr,Dr);
            [minval,minind] = min(abs(evals));
            [Y,inds] = sort(V(:,minind));
        end;
        Isort = inInds(inds);
        Wsort = W(inds,inds);
        Asort = A(inds,inds);
        clear A W D S Nr Dr V evals;
        [crit,affval,cutpoint] = choosecut(Wsort,Asort);
        if affval == 0
            carryon = (crit<thresh2);
        else
            carryon = (crit<thresh2)|((crit./affval)<thresh1);
        end
        
        finished = ~carryon;
        if(~finished);
            
            aff1 = Wsort(1:cutpoint,1:cutpoint);
            adj1 = Asort(1:cutpoint,1:cutpoint);
            aff2 = Wsort(cutpoint+1:end,cutpoint+1:end);
            adj2 = Asort(cutpoint+1:end,cutpoint+1:end);
            inds1 = Isort(1:cutpoint);
            inds2 = Isort(cutpoint+1:end);
            clear Asort Wsort Isort;
            
            queueitem{1} = aff1;queueitem{2} = adj1;queueitem{3} = inds1;  queue = {queue{:} queueitem};
            queueitem{1} = aff2;queueitem{2} = adj2;queueitem{3} = inds2;  queue = {queue{:} queueitem};
            
        else
            i1 = inInds(1);
            globalInds(inInds)=i1;
        end;
        
    elseif (numreg==2);  %best cut is trivial if only two regions
        if (W(2,1)>thresh2);   %slightly arbitrary choice of when to merge two single regions!!
            i1 = inInds(1);
            globalInds(inInds)=i1;
        end;
        
    end;
    
end
return;

function [breakval,mval,index] = choosecut(affreorder,adjace)
% choose based on mean between
N=size(affreorder,1)-1;
mmaffvals = zeros(N,1);

%Sparse version
sw = size(affreorder);
for yw = 1: sw(1)
    affreorder(yw,yw) = 0;
end


for k=1:N

    A1 = full(sum(sum(affreorder(1:k,1:k))));
    A2 = full(sum(sum(affreorder(k+1:end,k+1:end))));
    C = full(sum(sum(affreorder(1:k,k+1:end))));
    N1 = full(sum(sum(adjace(1:k,1:k))));
    N2 = full(sum(sum(adjace(k+1:end,k+1:end))));
    NC = full(sum(sum(adjace(1:k,k+1:end))));
    if A1>0
        A1= A1./N1;
    end
    if A2>0
        A2 = A2./N2;
    end
    if C>0
        C = C./NC;
    end
    mmaff = max(A1,A2);
    Cvals(k) = C;
    mmaffvals(k) = mmaff;
    
end

[breakval,index] = min(Cvals);
mval = mmaffvals(index);
return;
