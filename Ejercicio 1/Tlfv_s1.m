function [lfv,vsigma2]=Tlfv_s1(th)
global vyt
T=length(vyt);
mu=th(1);
delta1=th(3);
alpha1=th(4);
gamma1=th(5);
k=th(2);
e20=0.4; %std(vyt)^2;
e0=std(e20);
ind=e0<0;
sigma20=0.4; %std(vyt)^2;
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
