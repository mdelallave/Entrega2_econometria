function [lfv,vsigma2]=elfv(thtr)
% el modelo: y(t)=mu+z(t)*sigma(t), z(t)~N(0,1), e(t)=z(t)*sigma(t)
% sigma2(t)=k+alpha*e(t)^2+delta*sigma2(t)
%
global vyt
T=length(vyt);
% parametros iniciales originales
mu=thtr(1);
alpha1=exp(thtr(3));
beta1=exp(thtr(4))/(1+exp(thtr(4)));
alpha0=thtr(2);
gamma=-exp(thtr(5));
a0=0; %std(vyt);%0.4; %std(vyt)^2;
sigma20=std(vyt)^2; % varianza de los rendimientos

f=zeros(T,1);
vsigma2=zeros(T,1);
% Rellenamos la función de verosimilitud y el vector de varianzas
% estimadas usando un egarch este caso.
for t=1:T
    % varianza de las rendimientos.
    lsigma21=alpha0+alpha1*(1/(sigma20^0.5))*(abs(a0)+gamma*a0)+beta1*log(sigma20);
    sigma21=exp(lsigma21);
    vsigma2(t)=sigma21;
    
    f(t)=(1/((2*pi*sigma21)^0.5))*exp(-0.5*((vyt(t)-mu)^2)/sigma21);
    %fº densidad normal
    
    %innovaciones de los rendmientos en cada t.
    a=(vyt(t)-mu); 
    a0=a;
    sigma20=sigma21;
end
F=log(f);
lfv=ones(1,T)*F;
lfv=-lfv;

