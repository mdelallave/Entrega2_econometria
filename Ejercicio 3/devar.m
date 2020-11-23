function [varvec vardecm]=devar(N,M,MMA,bigV)
varvec=zeros(N*M,1);
vardecm=zeros(N*M,N);
vsigmak=zeros(N,N);
bigVx=bigV.^2;
for k=1:M
    vsigmak=vsigmak+MMA((k-1)*N+1:k*N,:)*bigVx*MMA((k-1)*N+1:k*N,:)';
    for i=1:N
        varvec((i-1)*M+k)=vsigmak(i,i);
    end
end
for j=1:N
    bigVxr=zeros(N,N);bigVxr(j,j)=bigVx(j,j);
    vsigmak=zeros(N,N);
    for k=1:M
        vsigmak=vsigmak+MMA((k-1)*N+1:k*N,:)*bigVxr*MMA((k-1)*N+1:k*N,:)';
        for i=1:N
            vardecm(k+(i-1)*M,j)=vsigmak(i,i)/varvec(k+(i-1)*M);
        end
    end
end
varvec=100*varvec;
vardecm=100*vardecm;