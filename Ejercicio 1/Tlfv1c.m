function [lfv,vsigma2]=Tlfv1c(thtr)
% el modelo: y(t)=mu+z(t)*sigma(t), z(t)~N(0,1), e(t)=z(t)*sigma(t) % z es
% epsilon y e es la a en los apuntes de clase
% sigma2(t)=k+alpha*e(t)^2+delta*sigma2(t)
%
global vyt
T=length(vyt);
mu=thtr(1);
k=exp(thtr(2)); % k es alpha_0 en los apuntes
delta1=exp(thtr(3))/(1+exp(thtr(3))+exp(thtr(4))+exp(thtr(5))); % delta es beta en los apuntes
alpha1=exp(thtr(4))/(1+exp(thtr(4))+exp(thtr(3))+exp(thtr(5)));
gamma1=2*exp(thtr(5))/(1+exp(thtr(4))+exp(thtr(3))+exp(thtr(5)));
e20=std(vyt)^2;%0.4; %std(vyt)^2;
e0=std(e20);
ind=e0<0;
sigma20=std(vyt)^2;%0.4; %std(vyt)^2;
f=zeros(T,1);
vsigma2=zeros(T,1);
for t=1:T
    sigma21=k+delta1*sigma20+(alpha1+gamma1*ind)*e20;vsigma2(t)=sigma21;
    f(t)=(1/((2*pi*sigma21)^0.5))*exp(-0.5*((vyt(t)-mu)^2)/sigma21);
    e2=(vyt(t)-mu)^2;e20=e2;
    sigma20=sigma21;
    e1=vyt(t)-mu; ind=e1<0;
end
F=log(f);
lfv=ones(1,T)*F;
lfv=-lfv;

