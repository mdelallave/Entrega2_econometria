function [lfv,vrho]=lfv_DCC(thtr)
%
%
global eps_est;
T=length(eps_est);
alpha=exp(thtr(1))/(1+exp(thtr(1))+exp(thtr(2)));
beta=exp(thtr(2))/(1+exp(thtr(1))+exp(thtr(2)));
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