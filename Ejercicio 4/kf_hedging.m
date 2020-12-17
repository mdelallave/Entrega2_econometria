function [lfv,vxit1t,vxitt,vxitT,MPtt,MPtT] = kf_hedging(th,ind)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global r1 r2
yt=r1;
T=length(yt);
xt=[ones(T,1) r2]';
%
betabar0=th(1);
betabar1=th(2);
if ind==1
    f11=th(3)/(1+abs(th(3)));
    f22=th(4)/(1+abs(th(4)));
    sigmaw=exp(-th(5))/10;
    sigma11=exp(-th(6))/10;
    sigma22=exp(-th(7))/10;
else
    f11=th(3);
    f22=th(4);
    sigmaw=th(5);
    sigma11=th(6);
    sigma22=th(7);
end
%Formando matrices
F=[f11 0;
    0  f22];
A=[betabar0;
   betabar1];
Q=[sigma11^2 0;
   0         sigma22^2];
R=sigmaw^2;
xit1t=zeros(2,1);
vxit1t=[];vxitt=[];MPtt=[];MPt1t=[];
vec_Pt1t=inv(eye(4)-kron(F,F))*Q(:);
Pt1t=reshape(vec_Pt1t,2,2);
lfv=log(2*pi)*T;
for t=1:T
    H=[1;
       r2(t)]; 
    wt=yt(t)-A'*xt(:,t)-H'*xit1t;
    Omegat=H'*Pt1t*H+R;
    Omegatinv=1/Omegat;
    lfv=lfv+0.5*log(Omegat)+0.5*(wt^2)/Omegat;
    xitt=xit1t+Pt1t*H*Omegatinv*wt;
    vxitt=[vxitt xitt];
    Ptt=Pt1t-Pt1t*H*Omegatinv*H'*Pt1t;
    MPtt=[MPtt Ptt];
    kt=F*Pt1t*H*Omegatinv;
    xit1t=F*xit1t+kt*wt;
    vxit1t=[vxit1t xit1t];
    Pt1t=(F-kt*H')*Pt1t*(F'-H*kt')+kt*R*kt'+Q;
    MPt1t=[MPt1t Pt1t];
end
%smoothing
vxitT=[];
MPtT=[];
vxiT=[vxitt(:,T) vxitT];
xitT=vxitt(:,T);
MPtT=[MPtt(:,2*(T-1)+1:2*(T-1)+2) MPtT];
PtT=MPtt(:,2*(T-1)+1:2*(T-1)+2);
for t=2:T
    Jt=MPtt(:,2*(T-t)+1:2*(T-t)+2)*F'*inv(MPt1t(:,2*(T-t)+1:2*(T-t)+2));
    xitT=vxitt(:,T-t+1)+Jt*(xitT-vxit1t(:,T-t+1));
    vxitT=[xitT vxitT];
    PtT=MPtt(:,2*(T-t)+1:2*(T-t)+2)*Jt*(PtT-MPt1t(:,2*(T-t)+1:2*(T-t)+2))*Jt';
    MPtT=[PtT MPtT];
end
%
end
