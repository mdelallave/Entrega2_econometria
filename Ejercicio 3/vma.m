function MMA=vma(MPHI,NAR,N,M)
ruido=zeros(N,M);
YT=zeros(N,M+NAR);
for i=1:N
    ruido(i,1)=1;
    for j=1:M
        YTt=ruido(:,j);
        for k=1:NAR
            YTt=YTt+MPHI(:,(k-1)*N+1:k*N)*YT(:,j+NAR-k);
        end
        YT(:,j+NAR)=YTt;
        MMAf(:,N*(j-1)+i)=YT(:,j+NAR);
    end
    ruido=zeros(N,M);
    YT=zeros(N,M+NAR);
end
MMA=[];  
for i=1:M
    MMA=[MMA;MMAf(:,(i-1)*N+1:i*N)];
end
    