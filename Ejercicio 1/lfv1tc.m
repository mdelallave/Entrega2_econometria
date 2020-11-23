function [lfv,vsigma2]=lfv1tc(thtr)
% el modelo: y(t)=mu+z(t)*sigma(t), z(t)~t(5) estandarizada, e(t)=z(t)*sigma(t)
% sigma2(t)=k+alpha*e(t)^2+delta*sigma2(t)
%
global vyt
T=length(vyt);
mu=thtr(1);
delta1=exp(thtr(3))/(1+exp(thtr(3))+exp(thtr(4)));
alpha1=exp(thtr(4))/(1+exp(thtr(4))+exp(thtr(3)));
k=exp(thtr(2));
v=exp(thtr(5));
e20=std(vyt)^2;%0.4; %std(vyt)^2;
sigma20=std(vyt)^2;%0.4; %std(vyt)^2;
f=zeros(T,1);
vsigma2=zeros(T,1);
for t=1:T
    sigma21=k+delta1*sigma20+alpha1*e20;vsigma2(t)=sigma21;
    f(t)=(gamma((v+1)/2)/(gamma(v/2)*(((v-2)*pi)^0.5)))*...
         (1/(sigma21^0.5))*((1+((vyt(t)-mu)^2)/((v-2)*sigma21))^(-(v+1)/2));
    e2=(vyt(t)-mu)^2;e20=e2;
    sigma20=sigma21;
end
F=log(f);
lfv=ones(1,T)*F;
lfv=-lfv;

