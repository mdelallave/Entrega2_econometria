function [lfv,vsigma2]=lfv2_garch11(thtr)
% el modelo: y(t)=mu+z(t)*sigma(t), z(t)~N(0,1), e(t)=z(t)*sigma(t)
% sigma2(t)=k+alpha*e(t)^2+delta*sigma2(t)
%
global r2;
T=length(r2);
mu=thtr(1);
delta1=exp(thtr(3))/(1+exp(thtr(3))+exp(thtr(4)));
alpha1=exp(thtr(4))/(1+exp(thtr(4))+exp(thtr(3)));
k=exp(thtr(2));
e20=std(r2)^2;%0.4; %std(r2)^2;
sigma20=std(r2)^2;%0.4; %std(r2)^2;
f=zeros(T,1);
vsigma2=zeros(T,1);
for t=1:T
    sigma21=k+delta1*sigma20+alpha1*e20;vsigma2(t)=sigma21;
    f(t)=(1/((2*pi*sigma21)^0.5))*exp(-0.5*((r2(t)-mu)^2)/sigma21);
    e2=(r2(t)-mu)^2;e20=e2;
    sigma20=sigma21;
end
F=log(f);
lfv=ones(1,T)*F;
lfv=-lfv;

