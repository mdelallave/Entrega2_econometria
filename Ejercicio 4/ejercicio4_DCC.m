clear;
global r1 r2 eps_est
datos1 = xlsread('data_WTI.xlsx','B4:C1721');
datos1=datos1(100:end,:);
r1=real(log(datos1(2:end,1)./datos1(1:end-1,1))); %rendimientos spot del oro
r2=real(log(datos1(2:end,2)./datos1(1:end-1,2))); %rendimientos futuro del oro 
r=[r1 r2];
T=length(r);
mu=mean(r);
vsigmat_1=(r1-mu(1)).^2;
vsigmat_2=(r2-mu(2)).^2;
vsigmat=[vsigmat_1 vsigmat_2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%regresión ratio de cobertura
X=[ones(T,1) r2];
betamco=inv(X'*X)*X'*r1;
rat_cob_mco=betamco(2);
Efectividad_mco=var(r1-rat_cob_mco*r2);
Efectividad_naive=var(r1-r2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paso 1: estimar un garch(1,1) para cada una de las dos series de
%rendimientos y obtener los residuos estandarizados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GARCH(1,1) para el spot
mu0=mu(1);
delta10=0.8; % beta
alpha10=0.1; % alpha_1
delta10tr=log(delta10/(1-delta10-alpha10)); % beta transformado
alpha10tr=log(alpha10/(1-alpha10-delta10)); % alpha_1 transformado
k0=(1-delta10-alpha10)*std(r1)^2; % alpha_0 (w)
k0tr=log(k0); % alpha_0 transformado
th0tr=[mu0,k0tr,delta10tr,alpha10tr]; % vector de param transformados
options=optimset('Display','off','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
thtropt=fminunc('lfv1_garch11',th0tr,options);% fminsearch  fminunc
mu1opt=thtropt(1);
delta1opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3)));
alpha1opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4)));
k1opt=exp(thtropt(2));
[lfv1,vsigma2_1]=lfv1_garch11(thtropt);
eps1=(r1-mu1opt)./(vsigma2_1.^0.5); % innovaciones estandarizadas para el precio spot
%
thopt1=[mu1opt k1opt delta1opt alpha1opt]'; 
x0=thopt1; n=length(x0);
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=diag(x0*1e-4);%auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfvs1_garch11',x0+auxi(:,i)+auxi(:,j))-...
                        feval('lfvs1_garch11',x0+auxi(:,i))-...
                        feval('lfvs1_garch11',x0+auxi(:,j))+...
                        feval('lfvs1_garch11',x0))/(auxi(j,j)^2);
        end
    end   
informd=inv(H0);
sgd1=sqrt(diag(informd));
%display('estimaciones de mu,alpha,delta y k y sus desviaciones típicas');
format long;
x1=x0;
%[x1 sgd1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GARCH(1,1) para el futuro
mu0=mu(2);
delta10=0.8;
alpha10=0.1;
delta10tr=log(delta10/(1-delta10-alpha10));
alpha10tr=log(alpha10/(1-alpha10-delta10));
k0=(1-delta10-alpha10)*std(r1)^2;
k0tr=log(k0);
th0tr=[mu0,k0tr,delta10tr,alpha10tr];
thtropt=fminunc('lfv2_garch11',th0tr,options);%fminsearch  fminunc
mu2opt=thtropt(1);
delta2opt=exp(thtropt(3))/(1+exp(thtropt(4))+exp(thtropt(3)));
alpha2opt=exp(thtropt(4))/(1+exp(thtropt(3))+exp(thtropt(4)));
k2opt=exp(thtropt(2));
[lfv2,vsigma2_2]=lfv2_garch11(thtropt);
eps2=(r2-mu2opt)./(vsigma2_2.^0.5);
%
thopt2=[mu2opt k2opt delta2opt alpha2opt]'; 
x0=thopt2; n=length(x0);
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=diag(x0*1e-4);%auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfvs2_garch11',x0+auxi(:,i)+auxi(:,j))-...
                        feval('lfvs2_garch11',x0+auxi(:,i))-...
                        feval('lfvs2_garch11',x0+auxi(:,j))+...
                        feval('lfvs2_garch11',x0))/(auxi(j,j)^2);
        end
    end   
informd=inv(H0);
sgd2=sqrt(diag(informd));
%display('estimaciones de mu,alpha,delta y k y sus desviaciones típicas');
x2=x0;
%[x2 sgd2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps_est=[eps1 eps2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Paso 2: estimar las correlaciones condicionales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha0=0.1;
beta0=0.75;
alpha0tr=log(alpha0/(1-alpha0-beta0));
beta0tr=log(beta0/(1-alpha0-beta0));
theta0tr=[alpha0tr beta0tr];
thetatropt=fminsearch('lfv_DCC',theta0tr,options);theta0tr=thetatropt; % Primera aproximación al óptimo
thetatropt=fminunc('lfv_DCC',theta0tr,options); % Segunda aproximación al óptimo a partir de la primera
alphaopt=exp(thetatropt(1))/(1+exp(thetatropt(1))+exp(thetatropt(2))); % Recupero el alpha original (deshago transformación)
betaopt=exp(thetatropt(2))/(1+exp(thetatropt(1))+exp(thetatropt(2))); % Recupero el beta original (deshago transformación)
[lfv_DCC,vrho]=lfv_DCC(thetatropt); % Calculamos las correlaciones óptimas cambiantes en el tiempo
cov12t=vrho.*(sqrt(vsigma2_1)).*(sqrt(vsigma2_2)); % Covarianzas óptimas cambiantes en el tiempo
thetaopt=[alphaopt betaopt]';
x0=thetaopt;n=length(x0);
 % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=diag(x0*1e-5);%auxi(2,2)=0.000000001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfvs_DCC',x0+auxi(:,i)+auxi(:,j))-...
                        feval('lfvs_DCC',x0+auxi(:,i))-...
                        feval('lfvs_DCC',x0+auxi(:,j))+...
                        feval('lfvs_DCC',x0))/(auxi(j,j)^2);
        end
    end   
informd=inv(H0);
sgd3=sqrt(diag(informd));
%display('estimaciones de mu,alpha,delta y k y sus desviaciones típicas');
x3=x0;
%[x3 sgd3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rat_cob_DCC = cov12t./vsigma2_2; % ratio de cobertura
Efectividad_DCC = var(r1-rat_cob_DCC.*r2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tablas
%Tabla Estimaciones
names_var={'mu_rend_spot    ';
           'alpha0_rend_spot';
           'beta_rend_spot  ';
           'alpha_rend_spot ';
           'mu_rend_fut     ';
           'alpha0_rend_fut ';
           'beta_rend_fut   ';
           'alpha_rend_fut  ';
           'alpha_DCC       ';
           'beta_DCC        '};
coefficient=[x1;x2;x3];
Std_Error=[sgd1;sgd2;sgd3];
t_statistic=coefficient./Std_Error;
p_value=1-normcdf(abs(t_statistic),0,1);p_value=2*p_value; %dos colas
disp('_________________________________________________________________________________________________________________')
disp('Coeficientes estimados')
TABLA1=table(coefficient,Std_Error,t_statistic,p_value,'RowNames',names_var)
%Tabla de efectividades del ratio de cobertura
names_var={'Efectividad naive';
           'Efectividad MCO  ';
           'Efectividad DCC  '};
coefficient=[Efectividad_naive;Efectividad_mco;Efectividad_DCC];
disp('_________________________________________________________________________________________________________________')
disp('_________________________________________________________________________________________________________________')
disp('Efectividades del ratio de cobertura')
TABLA2=table(coefficient,'RowNames',names_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gráficos
figure;
subplot(2,2,1)
plot(datos1(:,1));
title('Precio spot del crudo WTI (28/05/2014 a 16/11/2020)')
subplot(2,2,2)
plot(datos1(:,2));
title('Precio del futuro del crudo WTI (28/05/2014 a 16/11/2020)')
subplot(2,2,3)
plot(r1);
title('Rendimiento del spot del crudo WTI')
ylim([-0.12 0.1]);
subplot(2,2,4)
plot(r2);
title('Rendimiento del futuro del crudo WTI')
ylim([-0.12 0.1]);
%
% figure; % Varianzas condicionales spot y futuro
% tiledlayout(3,2)
% nexttile
% plot([vsigmat(:,1) vsigma2_1])
% title('Varianza condicional del rendimiento spot','FontSize',15);
% legend('r^2','varianza condicional estimada','FontSize',12);
% ylim([0 0.012]);
% 
% nexttile
% plot([vsigmat(:,2) vsigma2_2])
% title('Varianza condicional del rendimiento futuro','FontSize',15);
% legend('r^2','varianza condicional estimada','FontSize',12);
% ylim([0 0.012]);
% 
% nexttile
% plot(vsigma2_1)
% title('Varianza condicional del rendimiento spot','FontSize',15);
% ylim([0 0.0025]);
% 
% nexttile
% plot(vsigma2_2)
% title('Varianza condicional del rendimiento futuro','FontSize',15);
% ylim([0 0.002]);
% 
% nexttile([1 2])
% plot(cov12t,'LineWidth', 1.2)
% title('Covarianza condicional de los rendimientos spot y futuro','FontSize',15);
% ylim([0 0.004])%0.001]);

figure; % Varianzas condicionales spot y futuro
subplot(3,2,1)
plot([vsigmat(:,1) vsigma2_1])
title('Varianza condicional del rendimiento spot','FontSize',15);
legend('r^2','varianza condicional estimada','FontSize',12);
ylim([0 0.012]);

subplot(3,2,2)
plot([vsigmat(:,2) vsigma2_2])
title('Varianza condicional del rendimiento futuro','FontSize',15);
legend('r^2','varianza condicional estimada','FontSize',12);
ylim([0 0.012]);

subplot(3,2,3)
plot(vsigma2_1)
title('Varianza condicional del rendimiento spot','FontSize',15);
ylim([0 0.0025]);

subplot(3,2,4)
plot(vsigma2_2)
title('Varianza condicional del rendimiento futuro','FontSize',15);
ylim([0 0.002]);

subplot(3,2,[5 6])
plot(cov12t,'LineWidth', 1.2)
title('Covarianza condicional de los rendimientos spot y futuro','FontSize',15);
ylim([0 0.004])%0.001]);


figure; % Covarianza condicional y ratio de cobertura
plot([ones(length(rat_cob_DCC),1), rat_cob_mco.*ones(length(rat_cob_DCC),1), ...
    rat_cob_DCC],'LineWidth', 1.2);
title('Ratios de Cobertura','FontSize',18);
legend('naive', 'MCO', 'DCC','FontSize',15);
%ylim([-0.5 1.2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
