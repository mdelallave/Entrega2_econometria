%HAY DOS PROCESOS DE MARKOV INDEPENDIENTES: UNO PARA LAS BETAS DE LA
%REGRESIÓN Y OTRO PARA LAS VARIANZAS
clear;
global r1 r2
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
%% Markov Switching
beta_00=betamco(1);
beta_01=betamco(1)+0.001;
beta_10=betamco(2);
beta_11=betamco(2)+0.001;
sigma1=(((r1-X*betamco)'*(r1-X*betamco))/(length(r1)-2))^0.5;
sigma2=(((r1-X*betamco)'*(r1-X*betamco))/(length(r1)-2))^0.5+0.000001;
p11_1=0.95;
p22_1=0.95;
p11_2=0.9;
p22_2=0.9;
%
P=[p11_1*p11_2          p11_1*(1-p22_2)      (1-p22_1)*p11_2      (1-p22_1)*(1-p22_2);
   p11_1*(1-p11_2)      p11_1*p22_2          (1-p22_1)*(1-p11_2)  (1-p22_1)*p22_2;
   (1-p11_1)*p11_2      (1-p11_1)*(1-p22_2)  p22_1*p11_2          p22_1*(1-p22_2);
   (1-p11_1)*(1-p11_2)  (1-p11_1)*p22_2      p22_1*(1-p11_2)      p22_1*p22_2];
   
%Valores iniciales transformados
sigma1tr=log(sigma1);
sigma2tr=log(sigma2);
p11_1tr=log(p11_1/(1-p11_1));
p22_1tr=log(p22_1/(1-p22_1));
p11_2tr=log(p11_2/(1-p11_2));
p22_2tr=log(p22_2/(1-p22_2));
%
thtr0=[beta_00 beta_01 beta_10 beta_11 sigma1tr sigma2tr p11_1tr p22_1tr p11_2tr p22_2tr];
options=optimset('Display','iter','MaxFunEvals',10000,...
                 'MaxIter',10000,'TolFun',0.00001);
thtropt=fminsearch('lfv_MS_2',thtr0,options,1); %fminunc fminsearch
beta00opt=thtropt(1);
beta01opt=thtropt(2);
beta10opt=thtropt(3);
beta11opt=thtropt(4);
sigma1opt=exp(thtropt(5));
sigma2opt=exp(thtropt(6));
p11_1opt=exp(thtropt(7))/(1+exp(thtropt(7)));
p22_1opt=exp(thtropt(8))/(1+exp(thtropt(8)));
p11_2opt=exp(thtropt(9))/(1+exp(thtropt(9)));
p22_2opt=exp(thtropt(10))/(1+exp(thtropt(10)));
%
P=[p11_1opt*p11_2opt          p11_1opt*(1-p22_2opt)      (1-p22_1opt)*p11_2opt      (1-p22_1opt)*(1-p22_2opt);
   p11_1opt*(1-p11_2opt)      p11_1opt*p22_2opt          (1-p22_1opt)*(1-p11_2opt)  (1-p22_1opt)*p22_2opt;
   (1-p11_1opt)*p11_2opt      (1-p11_1opt)*(1-p22_2opt)  p22_1opt*p11_2opt          p22_1opt*(1-p22_2opt);
   (1-p11_1opt)*(1-p11_2opt)  (1-p11_1opt)*p22_2opt      p22_1opt*(1-p11_2opt)      p22_1opt*p22_2opt];
thopt=[beta00opt beta01opt beta10opt beta11opt sigma1opt sigma2opt ...
       p11_1opt p22_1opt p11_2opt p22_2opt]';
x0=thopt; n=length(x0);
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=diag(x0*1e-4);auxi(2,2)=0.0001;%auxi(7,7)=0.00001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('lfv_MS_2',x0+auxi(:,i)+auxi(:,j),0)-...
                        feval('lfv_MS_2',x0+auxi(:,i),0)-...
                        feval('lfv_MS_2',x0+auxi(:,j),0)+...
                        feval('lfv_MS_2',x0,0))/(auxi(j,j)^2);
        end
    end   
informd=inv(H0);
sgd=sqrt((diag(informd)));
%display('estimaciones de mu,alpha,delta y k y sus desviaciones típicas');
format long;
%[x0 sgd]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lfv,xitt,xit1t]=lfv_MS_2(x0,0);
%
%smoothing (Kim's algorithm)
xitT=zeros(T-1,4);
xitT=[xitT; xitt(T,:)];
for t=1:T-1
    a=xitt(T-t,:)'.*(P'*(xitT(T-t+1,:)'./xit1t(T+1-t,:)'));
    xitT(T-t,:)=a';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio_cobtt=(xitt(:,1)+xitt(:,2)).*beta10opt+(xitt(:,3)+xitt(:,4)).*beta11opt;
ratio_cobtT=(xitT(:,1)+xitT(:,2)).*beta10opt+(xitT(:,3)+xitT(:,4)).*beta11opt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Efectividad_MS_tt=var(r1-ratio_cobtt.*r2);
Efectividad_MS_tT=var(r1-ratio_cobtT.*r2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tabla de estimaciones
names_var={'beta0-régimen 1 ';
           'beta0-régimen 2 ';           
           'beta1-régimen 1 ';
           'beta1-régimen 2 ';           
           'sigma-régimen 1 ';
           'sigma-régimen 2 ';
           'p11_1           ';
           'p22_1           ';
           'p11_2           ';
           'p22_2           '};
coefficient=x0;
Std_Error=sgd;
t_statistic=coefficient./Std_Error;
p_value=1-normcdf(abs(t_statistic),0,1);p_value=2*p_value; %dos colas
disp('_________________________________________________________________________________________________________________')
disp('Coeficientes estimados')
TABLA1=table(coefficient,Std_Error,t_statistic,p_value,'RowNames',names_var)
%Tabla de efectividades del ratio de cobertura
names_var={'Efectividad naive  ';
           'Efectividad MCO    ';
           'Efectividad MS_{t|T}';
           'Efectividad MS_{t|t}'};
coefficient=[Efectividad_naive;Efectividad_mco;Efectividad_MS_tT;Efectividad_MS_tt];
disp('_________________________________________________________________________________________________________________')
disp('_________________________________________________________________________________________________________________')
disp('Efectividades del ratio de cobertura')
TABLA2=table(coefficient,'RowNames',names_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(3,2,1);
plot([xitt(:,1) xitT(:,1)]);
title('Probabilidad de que una observación esté gobernada por el régimen 1');
legend('xi_{t|t}','xi_{t|T}');
subplot(3,2,2);
plot([xitt(:,2) xitT(:,2)]);
title('Probabilidad de que una observación esté gobernada por el régimen 2');
legend('xi_{t|t}','xi_{t|T}');
subplot(3,2,3);
plot([xitt(:,3) xitT(:,3)]);
title('Probabilidad de que una observación esté gobernada por el régimen 3');
legend('xi_{t|t}','xi_{t|T}');
subplot(3,2,4);
plot([xitt(:,4) xitT(:,4)]);
title('Probabilidad de que una observación esté gobernada por el régimen 4');
legend('xi_{t|t}','xi_{t|T}');
subplot(3,2,5);
plot([ones(length(ratio_cobtt),1), rat_cob_mco.*ones(length(ratio_cobtt),1)...
      ratio_cobtt ratio_cobtT]);
title('Ratios de Cobertura');
legend('naive', 'MCO', 'ratio-cob_{t|t}', 'ratio-cob_{t|T}');