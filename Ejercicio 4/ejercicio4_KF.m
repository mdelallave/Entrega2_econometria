clear;clc;
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
%Kalman Filter
% ESTOS PARÁMETROS INICIALES NO SON NECESARIOS SI USAMOS LOS VALORES
% INICIALES DEL PROFESOR. BORRAR EN CASO DE NO USARLOS.
betabar00=betamco(1);
betabar10=betamco(2);
f110=0.5;f110tr=log(f110/(1-f110));
f220=0.5;f220tr=log(f220/(1-f220));
sigmaw0=std(r1-X*betamco);sigmaw0tr=log(sigmaw0);
sigma110=0.001;sigma110tr=log(sigma110);
sigma220=0.1;sigma220tr=log(sigma220);
%
%th0tr=[betabar00 betabar10 2.3 2.3 2 3 3]; NO FUNCIONA
th0tr=[0.000071 0.8665 -0.8 0.93/(1-0.93) -log(10*0.000058) -log(0.066) -log(0.335)]; % VALORES INICIALES PROFESOR
%th0tr=[-betabar00 betabar10 f110tr f220tr -sigmaw0tr -sigma110tr -sigma220tr]; % VALORES INICIALES PROPUESTOS (MANU)

%
options=optimset('Display','Iter','MaxFunEvals',10000,'MaxIter',10000,... 
                 'TolFun',0.00001);
%             
[lfv,vxit1t,vxitt,vxitT,MPtt,MPtT] = kf_hedging(th0tr,1);
%
thtropt=fminunc('kf_hedging',th0tr,options,1); %fminsearch  fminunc
th0tr=thtropt; thtropt=fminunc('kf_hedging',th0tr,options,1);
%
betabar0=thtropt(1);
betabar1=thtropt(2);
f11=thtropt(3)/(1+abs(thtropt(3)));
f22=thtropt(4)/(1+abs(thtropt(4)));
sigmaw=exp(-thtropt(5))/10;
sigma11=exp(-thtropt(6))/10;
sigma22=exp(-thtropt(7))/10;
%
thopt=[betabar0 betabar1 f11 f22 sigmaw sigma11 sigma22]';
x0=thopt; n=length(x0);
    % Calculando el hessiano numéricamente
    H0=zeros(n,n);
    auxi=diag(x0*1e-4);%auxi(4,4)=0.0001;auxi(7,7)=0.00001;
    for i=1:n
        for j=1:n
            H0(i,j)=(feval('kf_hedging',x0+auxi(:,i)+auxi(:,j),0)-...
                        feval('kf_hedging',x0+auxi(:,i),0)-...
                        feval('kf_hedging',x0+auxi(:,j),0)+...
                        feval('kf_hedging',x0,0))/(auxi(j,j)^2);
        end
    end   
informd=inv(H0);
sgd1=sqrt(abs(diag(informd)));
%display('estimaciones de mu,alpha,delta y k y sus desviaciones típicas');
format long;
[x0 sgd1]
%
[lfv,xit1t,xitt,xitT,MPtt,MPtT] = kf_hedging(thopt,0);
ratio_cobtt=xitt(2,2:end)'+betabar1;
ratio_cobtT=xitT(2,:)'+betabar1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot([ones(length(ratio_cobtt),1), rat_cob_mco.*ones(length(ratio_cobtt),1)...
      ratio_cobtt ratio_cobtT]);
title('Ratios de Cobertura');
legend('naive', 'MCO', 'ratio-cob_{t|t}', 'ratio-cob_{t|T}');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Efectividad_KF_tT=var(r1(2:end)-ratio_cobtT.*r2(2:end));
Efectividad_KF_tt=var(r1(2:end)-ratio_cobtt.*r2(2:end));
%
%Tabla Estimaciones
names_var={'beta_bar0 ';
           'beta_bar1 ';           
           'f11       ';
           'f22       ';
           'sigmaw    ';
           'sigma11   ';
           'sigma22   '};
coefficient=x0;
Std_Error=sgd1;
t_statistic=coefficient./Std_Error;
p_value=1-normcdf(abs(t_statistic),0,1);p_value=2*p_value; %dos colas
disp('_________________________________________________________________________________________________________________')
disp('Coeficientes estimados')
TABLA1=table(coefficient,Std_Error,t_statistic,p_value,'RowNames',names_var)
%Tabla de efectividades del ratio de cobertura
names_var={'Efectividad naive  ';
           'Efectividad MCO    ';
           'Efectividad KF_{t|T}';
           'Efectividad KF_{t|t}'};
coefficient=[Efectividad_naive;Efectividad_mco;Efectividad_KF_tT;Efectividad_KF_tt];
disp('_________________________________________________________________________________________________________________')
disp('_________________________________________________________________________________________________________________')
disp('Efectividades del ratio de cobertura')
TABLA2=table(coefficient,'RowNames',names_var)

