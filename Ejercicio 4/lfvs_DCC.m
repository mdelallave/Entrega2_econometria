function [lfv,vrho]=lfvs_DCC(th)
%
%
global eps_est;
T=length(eps_est);
alpha=th(1);
beta=th(2);
Qbar=corr(eps_est);
Qt=Qbar;
vrho=zeros(T,1);
lfv=0;
for t=1:T
    Qt=(1-alpha-beta).*Qbar+alpha.*eps_est(t,:)'*eps_est(t,:)+...
       beta.*Qt;
    Qtilde=diag(sqrt(diag(Qt)));
    Rt=(Qtilde\Qt)/Qtilde;
    vrho(t)=Rt(1,2);
    lfv=lfv+log(det(Rt))+(eps_est(t,:)/Rt)*eps_est(t,:)';
end
lfv=0.5*lfv;