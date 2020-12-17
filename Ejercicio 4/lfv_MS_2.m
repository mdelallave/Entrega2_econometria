function [lfv,xitt,xit1t]=lfv_MS_2(thtr,ind)
% parametros transformados
global r1 r2
T=length(r1);
Z=[ones(length(r1),1) r2];
thtr=real(thtr);
N=4;
%
beta00=thtr(1);
beta01=thtr(2);
beta10=thtr(3);
beta11=thtr(4);
if ind==1
    sigma1=exp(thtr(5));
    sigma2=exp(thtr(6));
    p11_1=exp(thtr(7))/(1+exp(thtr(7)));
    p22_1=exp(thtr(8))/(1+exp(thtr(8)));
    p11_2=exp(thtr(9))/(1+exp(thtr(9)));
    p22_2=exp(thtr(10))/(1+exp(thtr(10)));
else
    sigma1=thtr(5);
    sigma2=thtr(6);
    p11_1=thtr(7);
    p22_1=thtr(8);
    p11_2=thtr(9);
    p22_2=thtr(10);
end
beta0=[beta00 beta10]';
beta1=[beta01 beta11]';
% computing eta
eta=zeros(T,N);
for t=1:T
    eta(t,1)=(1/(sigma1*(2*pi)^0.5))*exp(-((r1(t)-Z(t,:)*beta0).^2)/(2*sigma1^2));
    eta(t,2)=(1/(sigma2*(2*pi)^0.5))*exp(-((r1(t)-Z(t,:)*beta0).^2)/(2*sigma2^2));
    eta(t,3)=(1/(sigma1*(2*pi)^0.5))*exp(-((r1(t)-Z(t,:)*beta1).^2)/(2*sigma1^2));
    eta(t,4)=(1/(sigma2*(2*pi)^0.5))*exp(-((r1(t)-Z(t,:)*beta1).^2)/(2*sigma2^2));
end
% computing xi(t|t) and xi(t+1|t)
    % computing ergodic probabilities
   P=[p11_1*p11_2          p11_1*(1-p22_2)      (1-p22_1)*p11_2      (1-p22_1)*(1-p22_2);
      p11_1*(1-p11_2)      p11_1*p22_2          (1-p22_1)*(1-p11_2)  (1-p22_1)*p22_2;
      (1-p11_1)*p11_2      (1-p11_1)*(1-p22_2)  p22_1*p11_2          p22_1*(1-p22_2);
      (1-p11_1)*(1-p11_2)  (1-p11_1)*p22_2      p22_1*(1-p11_2)      p22_1*p22_2];
    A=[eye(N)-P; ones(1,N)];
    eN1=[zeros(N,1);
        1];
    ppi=inv(A'*A)*A'*eN1;
%    
xi10=ppi';
xitt=zeros(T,N);
xit1t=zeros(T,N);
for t=1:T
    xitt(t,:)=(xi10.*eta(t,:))/((xi10.*eta(t,:))*ones(N,1));
    xit1t(t,:)=xitt(t,:)*P';
    xi10=xit1t(t,:);
end
xit1t=[xi10;xit1t];
f=(xit1t(1:T,:).*eta)*ones(N,1);
lf=log(f);
lfv=ones(1,T)*lf;
%
lfv=-lfv;
