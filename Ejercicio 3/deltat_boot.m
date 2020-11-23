function MDeltat=deltat_boot(resid,yhat,nboot)
global Y NAR M H N vnames
%
MDeltat=zeros(N,N,nboot);
for jj=1:nboot
    jj
    nobs=length(resid);
    ystar=yhat+resid(ceil(size(resid,1)*rand(nobs,1)),:);
    %size(ystar)
    ystar=[Y(1:NAR,:);ystar];
    results=fvar(ystar,NAR);
    MPHI=mphi(results,NAR,N);
    %Calculando la matriz de varianzas y covarianzas
    for i = 1:N
        err(:,i) = results(i).resid;
    end
    Omega = err'*err/rows(err);
    P=eye(N);
    [MFRI,MMA]=fri(results,MPHI,NAR,N,M,P,vnames);
    DeltaC=zeros(N,N);
    for i=1:N
        for j=1:N
            num=0;
            den=0;
            for h=1:H
                ei=zeros(N,1);ei(i)=1;
                ej=zeros(N,1);ej(j)=1;
                num=num+(ei'*MMA(1+N*(h-1):N*h,:)*Omega*ej)^2;
                den=den+ei'*MMA(1+N*(h-1):N*h,:)*Omega*MMA(1+N*(h-1):N*h,:)'*ei;
            end
            DeltaC(i,j)=(1/Omega(j,j))*num/den;
        end
    end
    Deltat=zeros(N,N);
    for i=1:N
        for j=1:N
            S=sum(DeltaC');
            Deltat(i,j)=DeltaC(i,j)/S(i);
        end
    end
    Deltat=Deltat.*100;
    MDeltat(:,:,jj)=Deltat;
end

