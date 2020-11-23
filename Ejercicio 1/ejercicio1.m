%% Limpiamos
close('all')
clearvars
clc

%% Generamos datos con una estructura determinada

rng(1585) % Semilla para replicar los datos

global vyt
T=500;
mu=0.01;
phi=0.96;
a=4;
b=0.2;
vsigma2t=[];
vyt=[];
sigma2t0=a*b;
for t=1:T
    u=rand(1,1);
    if u<phi
        sigma2t=sigma2t0;
        yt=mu+sigma2t*randn(1);
        vsigma2t=[vsigma2t;sigma2t];
        vyt=[vyt;yt];
        sigma2t0=sigma2t;
    else
        sigma2t=gamrnd(a,b);
        yt=mu+sigma2t*randn(1);
        vsigma2t=[vsigma2t;sigma2t];
        vyt=[vyt;yt];
        sigma2t0=sigma2t;
    end
end

figure(1);
subplot(2,1,1);
plot(vyt);
title('niveles');
subplot(2,1,2);
plot(vsigma2t);
title('varianza');

%% Construimos un modelo GARCH(1,1) con innovaciones normales (g11)
% Valores iniciales 
mu0=mu;
delta10=0.8;
alpha10=0.12;

% Transformamos los parámetros
delta10tr=log(delta10/(1-delta10-alpha10));
alpha10tr=log(alpha10/(1-alpha10-delta10));
k0=(1-delta10-alpha10)*0.4;
k0tr=log(k0);

% Optimizamos los parametros transformados
th0tr=[mu0,k0tr,delta10tr,alpha10tr];
options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
thtropt=fminsearch('lfv1c',th0tr,options); 

% Recuperamos los parámetros, ahora óptimos
muopt=thtropt(1);
delta1opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3)));
alpha1opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4)));
kopt=exp(thtropt(2));

% Obtenemos la volatilidad condicional y maximizamos la función de 
% verosimilitud para los parametros óptimos.
[lfv_g11,vsigma2_g11]=lfv1c(thtropt);
figure(2);
plot([vsigma2t vsigma2_g11]);

thopt_g11=[muopt kopt delta1opt alpha1opt]'; % Vector de parametros optimos
x0_g11=thopt_g11; 
n=length(x0_g11);

    % Calculando el gradiente numéricamente
    auxi=1e-7*eye(n);
    Grad0=zeros(n,1);
    for i=1:n
        Grad0(i)=(feval('lfv_s1',x0_g11+auxi(:,i)) - ...
            feval('lfv_s1',x0_g11-auxi(:,i)))/(2*auxi(i,i));
    end
    
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=1e-5*eye(n);
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfv_s1',x0_g11+auxi(:,i)+auxi(:,j))-...
                        feval('lfv_s1',x0_g11+auxi(:,i))-...
                        feval('lfv_s1',x0_g11+auxi(:,j))+...
                        feval('lfv_s1',x0_g11))/(auxi(j,j)^2);
        end
    end
 
    
informd_g11=inv(H0);
sgd_g11=sqrt(diag(informd_g11));

[x0_g11 sgd_g11]

% Obtenemos medida de calidad del ajuste
[aic_g11,bic_g11] = aicbic(lfv_g11(end),length(thopt_g11),T);

%% EGARCH(1,1) con innovación normal (eg11)
% Valores iniciales
mu0=mu;
alpha00=-0.5; 
alpha10=0.2;
beta10=0.9;
gamma0=-0.2;

% Valores transformados
beta10tr=log(beta10/(1-beta10)); 
alpha10tr=log(alpha10);
gamma0tr=log(-gamma0);
th0tr=[mu0,alpha00,alpha10tr,beta10tr,gamma0tr];

% Obtenemos los parámetros optimos usando un egarch y fº de
% versomilitud basada en la norma
options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
thtropt=fminunc('elfv1c',th0tr,options);

muopt=thtropt(1);
alpha1opt=exp(thtropt(3));
beta1opt=exp(thtropt(4))/(1+exp(thtropt(4)));
alpha0opt=thtropt(2);
gammaopt=-exp(thtropt(5));

% Con los parámetros óptimos, calculamos volatilidad condicional y valor de
% fº verosimilitud
[lfv_eg11,vsigma2_eg11]=elfv1c(thtropt);
figure(3);
plot([vsigma2t vsigma2_eg11]);

thopt_eg11=[muopt alpha0opt alpha1opt beta1opt gammaopt]'; 
x0_eg11=thopt_eg11; 
n=length(x0_eg11);

 % Calculando el gradiente numéricamente
    auxi=1e-7*eye(n); auxi(2,2)=0.000000001;
    Grad0=zeros(n,1);
    for i=1:n
        Grad0(i)=(feval('elfv_s1',x0_eg11+auxi(:,i))-...
            feval('elfv_s1',x0_eg11-auxi(:,i)))/(2*auxi(i,i));
    end
    
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=1e-5*eye(n); %auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('elfv_s1',x0_eg11+auxi(:,i)+auxi(:,j))-...
                        feval('elfv_s1',x0_eg11+auxi(:,i))-...
                        feval('elfv_s1',x0_eg11+auxi(:,j))+...
                        feval('elfv_s1',x0_eg11))/(auxi(j,j)^2);
        end
    end
 
    
informd_eg11=inv(H0);
sgd_eg11=sqrt(diag(informd_eg11));

[x0_eg11 sgd_eg11]

% Obtenemos medida de calidad del ajuste
[aic_eg11,bic_eg11] = aicbic(lfv_eg11(end),length(thopt_eg11),T);

%% Garch(1,1) en media con innovaciones normales (mg11)
% Condiciones iniciales
mu0=mu;
c0=0; % Condicion inicial: no prima de riesgo.
delta10=0.5;
alpha10=0.3;

% Transformacion condiciones iniciales
delta10tr=log(delta10/(1-delta10-alpha10));
alpha10tr=log(alpha10/(1-alpha10-delta10));
k0=(1-delta10-alpha10)*std(vyt)^2;
k0tr=log(k0);
th0tr=[mu0,k0tr,delta10tr,alpha10tr,c0]

% Optimizamos los parámetros
options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
thtropt=fminunc('lfv1mc',th0tr,options);
muopt=thtropt(1);
delta1opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3)));
alpha1opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4)));
kopt=exp(thtropt(2));
copt=thtropt(5);

% Calculamos vol condicional y valor de fº verosimilitud
[lfv_mg11,vsigma2_mg11]=lfv1mc(thtropt);
figure(4);
plot([vsigma2t vsigma2_mg11]);


thopt_mg11=[muopt kopt delta1opt alpha1opt copt]'; 
x0_mg11=thopt_mg11; n=length(x0_mg11);
    % Calculando el gradiente numéricamente
    auxi=1e-7*eye(n); auxi(2,2)=0.000000001;
    Grad0=zeros(n,1);
    for i=1:n
        Grad0(i)=(feval('lfvm_s1',x0_mg11 + auxi(:,i))-...
            feval('lfvm_s1',x0_mg11 - auxi(:,i)))/(2*auxi(i,i));
    end
    
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=1e-5*eye(n);auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfvm_s1',x0_mg11+auxi(:,i)+auxi(:,j))-...
                        feval('lfvm_s1',x0_mg11+auxi(:,i))-...
                        feval('lfvm_s1',x0_mg11+auxi(:,j))+...
                        feval('lfvm_s1',x0_mg11))/(auxi(j,j)^2);
        end
    end
 
    
informd_mg11=inv(H0);
sgd_mg11=sqrt(diag(informd_mg11));

[x0_mg11 sgd_mg11]

[aic_mg11,bic_mg11] = aicbic(lfv_mg11(end),length(thopt_mg11),T);

%% GARCH(1,1) con innovaciones t-student (tg11)
% Parámetros iniciales
mu0=mu;
delta10=0.5;
alpha10=0.3;
gradosl0=5; % Condicion inicial para grados de libertad

% Transformamos los parametros iniciales
delta10tr=log(delta10/(1-delta10-alpha10));
alpha10tr=log(alpha10/(1-alpha10-delta10));
k0=(1-delta10-alpha10)*std(vyt)^2;
k0tr=log(k0);
gradosl0tr=log(gradosl0);

% Escribimos el vector de parametros transformados y lo optimizamos
th0tr=[mu0,k0tr,delta10tr,alpha10tr,gradosl0tr];
options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);

% Minimizamos la función de una tstudent
thtropt=fminunc('lfv1tc',th0tr,options);

% Recuperamos los valores en niveles, pero ahora optimos.
muopt=thtropt(1);
delta1opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3)));
alpha1opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4)));
kopt=exp(thtropt(2));
gradoslopt=exp(thtropt(5));
[lfv_tg11,vsigma2_tg11]=lfv1tc(thtropt);
figure(5);
plot([vsigma2t vsigma2_tg11]);

% Calculamos la matriz de varianzas-covarianzas
thopt_tg11=[muopt kopt delta1opt alpha1opt gradoslopt]'; 
x0_tg11=thopt_tg11; n=length(x0_tg11);

    % Calculando el gradiente numéricamente
    auxi=1e-7*eye(n); auxi(2,2)=0.000000001;
    Grad0=zeros(n,1);
    for i=1:n
        Grad0(i)=(feval('lfvt_s1',x0_tg11+auxi(:,i))-...
            feval('lfvt_s1',x0_tg11-auxi(:,i)))/(2*auxi(i,i));
    end
    
    % Calculando el hessiano numéricamente
    H0 = zeros(n,n);
    auxi=1e-5*eye(n);auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfvt_s1',x0_tg11+auxi(:,i)+auxi(:,j))-...
                        feval('lfvt_s1',x0_tg11+auxi(:,i))-...
                        feval('lfvt_s1',x0_tg11+auxi(:,j))+...
                        feval('lfvt_s1',x0_tg11))/(auxi(j,j)^2);
        end
    end 
    
informd_tg11=inv(H0);
sgd_tg11=sqrt(diag(informd_tg11));

[x0_tg11 sgd_tg11]

[aic_tg11,bic_tg11] = aicbic(lfv_tg11(end),length(thopt_tg11),T);

%% TGARCH(1,1) (nombramos a las variables usando tga11)
% Condiciones iniciales
mu0=mu;
delta10=0.5;
alpha10=0.3;
gamma10=0.1;

% Transformamos parametros
delta10tr=log(delta10/(1-delta10-alpha10-gamma10/2));
alpha10tr=log(alpha10/(1-alpha10-delta10-gamma10/2));
gamma10tr=log(0.5*gamma10/(1-alpha10-delta10-gamma10/2));
k0=(1-delta10-alpha10-gamma10/2)*std(vyt)^2;
k0tr=log(k0);
th0tr=[mu0,k0tr,delta10tr,alpha10tr,gamma10tr];

% Optimizamos los parámetros minimizando la función de verosimilitud
options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
thtropt=fminunc('Tlfv1c',th0tr,options);

muopt=thtropt(1);
delta1opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3))+exp(thtropt(5)));
alpha1opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4))+exp(thtropt(5)));
gamma1opt=2*exp(thtropt(5))/(1+exp(thtropt(3))+exp(thtropt(4))+exp(thtropt(5)));
kopt=exp(thtropt(2));

% Optimizamos fº de verosimilitud y obtenemos vol condicional con
% parametros óptimos
[lfv_tga11,vsigma2_tga11]=Tlfv1c(thtropt);
figure(6);
plot([vsigma2t vsigma2_tga11]);

thopt_tga11=[muopt kopt delta1opt alpha1opt gamma1opt]'; 
x0_tga11=thopt_tga11; n=length(x0_tga11);

    % Calculando el gradiente numéricamente
    auxi=1e-7*eye(n); auxi(2,2)=0.000000001;
    Grad0=zeros(n,1);
    for i=1:n
        Grad0(i)=(feval('Tlfv_s1',x0_tga11+auxi(:,i))-...
            feval('Tlfv_s1',x0_tga11-auxi(:,i)))/(2*auxi(i,i));
    end
    
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=1e-5*eye(n); auxi(2,2)=0.000000001;
    auxi(4,4)=(1/1000)*auxi(4,4); % Evitamos num negativos en la varianza
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('Tlfv_s1',x0_tga11+auxi(:,i)+auxi(:,j))-...
                        feval('Tlfv_s1',x0_tga11+auxi(:,i))-...
                        feval('Tlfv_s1',x0_tga11+auxi(:,j))+...
                        feval('Tlfv_s1',x0_tga11))/(auxi(j,j)^2);
        end
    end
 
    
informd_tga11=inv(H0);
sgd_tga11=sqrt(diag(informd_tga11));

[x0_tga11 sgd_tga11]

[aic_tga11,bic_tga11] = aicbic(lfv_tga11(end),length(thopt_tga11),T);

%% Comparamos los modelos
AIC = [aic_g11 aic_eg11 aic_mg11 aic_tg11 aic_tga11]';
BIC = [bic_g11 bic_eg11 bic_mg11 bic_tg11 bic_tga11]';

names_var = {'GARCH ';
            'EGARCH';
            'GARCHM';
            'GARCHt';
            'TGARCH'};
       
tabla1 = table(AIC, BIC,'RowNames',names_var)