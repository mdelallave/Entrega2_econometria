function MPHI=mphi(results,NAR,N)
A=[];
for i=1:N
    A=[A;results(i).beta(1:N*NAR)'];
end

MPHI=zeros(N,N*NAR);
for i=1:NAR
    PHI=zeros(N,N);
    for j=1:N
        PHI(:,j)=A(:,(j-1)*NAR+i);
    end
    MPHI(:,(i-1)*N+1:i*N)=PHI;
end
    