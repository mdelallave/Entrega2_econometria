function [lfv,vsigma2]=lfvt_s1(th)
global vyt
T=length(vyt);
mu=th(1);
delta1=th(3);
alpha1=th(4);
k=th(2);
v=th(5);
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
