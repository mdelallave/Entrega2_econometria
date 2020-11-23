%% Ejercicio 3
clear;
close all
clc

global vyt Y T NAR M H N vnames
%% Cargamos los datos %%
data = xlsread('data','data', 'B2:G1045');

% Sacamos rendimientos logaritmicos para cada empresa
r(:,:) = log(data((2:end),:) ./ data((1:end-1),:));
mu = mean(r);
vsigma2t = (r - mu).^2;

%% AJUSTAMOS UN EGARCH con innovaciones normales (Recoge efecto leverage)
% Fijamos valores iniciales genericos

alpha00 = -0.5; 
alpha10 = 0.2;
beta10 = 0.9;
gamma0 = -0.2;
beta10tr = log(beta10/(1-beta10)); 
alpha10tr = log(alpha10);
gamma0tr = log(- gamma0);
T=length(r);

% Creamos una matriz para los elementos del bucle
paramopt = zeros(5,5);
vsigma2 = zeros(T,5);
lfv = zeros(T,5);

for i=1:size(r,2)
    vyt = r(:,i);
    mu0 = mu(:,i);
    th0tr = [mu0, alpha00, alpha10tr ,beta10tr, gamma0tr];
    options=optimset('Display','iter','MaxFunEvals',10000,...
    'MaxIter',10000,'TolFun',0.00001);
    thtropt = fminunc('elfv',th0tr,options);

    paramopt(i,1) = thtropt(1);
    paramopt(i,2) = exp(thtropt(3));
    paramopt(i,3) = exp(thtropt(4))/(1+exp(thtropt(4)));
    paramopt(i,4) = thtropt(2);
    paramopt(i,5) = -exp(thtropt(5));

    % Con los parámetros óptimos, calculamos volatilidad condicional y 
    % valor de fº verosimilitud
    [lfv(:,i), vsigma2(:,i)] = elfv(thtropt);
end 
vsigma2_ib35 = vsigma2(:,1); vsigma2_sab = vsigma2(:,2);
vsigma2_bbva = vsigma2(:,3); vsigma2_san = vsigma2(:,4);
vsigma2_IB = vsigma2(:,5); vsigma2_bka = vsigma2(:,6);

figure;
subplot(3,2,1)
plot([vsigma2t(:,1) vsigma2(:,1)]);
    legend('Varianza serie Ibex 35', 'Ajuste EGARCH(1,1) Ibex 35')
subplot(3,2,2)
plot([vsigma2t(:,2) vsigma2(:,2)]);
    legend('Varianza serie Sabadell', 'Ajuste EGARCH(1,1) Sabadell')
subplot(3,2,3)
plot([vsigma2t(:,3) vsigma2(:,3)]);
    legend('Varianza serie BBVA', 'Ajuste EGARCH(1,1) BBVA')
subplot(3,2,4)
plot([vsigma2t(:,4) vsigma2(:,4)]);
    legend('Varianza serie Santander', 'Ajuste EGARCH(1,1) Santander')
subplot(3,2,5)
plot([vsigma2t(:,5) vsigma2(:,5)]);
    legend('Varianza serie Índice Bancario', 'Ajuste EGARCH(1,1) Índice Bancario')
subplot(3,2,6)
plot([vsigma2t(:,6) vsigma2(:,6)]);
    legend('Varianza serie Bankia', 'Ajuste EGARCH(1,1) Bankia')
    
    
%% Estimamos un VAR a los logaritmos de las volatilidades estimadas.

Y = [log(vsigma2)];

% Contratamos un test sobre razón de verosimilitud para determinar el orden
% óptimo del VAR.
sims = 1;
maxlag = 16;
minlag = 1;
lrratio(Y,maxlag,minlag,sims);

% Contrastamos si elegimos un retardo (p0) o 14 retardos (p1)
% Contrastar p1 retardos contra p0
display('_______________________________________________________________');
display('Contraste 14 retardos versus 1');
likratio(Y,14,1,0.05,1)
display('_______________________________________________________________');

% Rechazamos H_0: nos quedamos con un VAR(1).

%%  Estimamos el modelo VAR %%
T = length(Y);  %tamaño muestral
NAR = 1;  % Número de retardos del VAR
M = 50;   % Número de retardos del VMA(inf)
H = 12;   % Horizonte de predicción para la descomposición de la varianza
N = size(Y,2);

[results]=fvar(Y,NAR);
MPHI=mphi(results,NAR,N);

% Matriz de varianzas y covarianzas
for i = 1:N
    err(:,i) = results(i).resid;
end
Omega = err'*err/rows(err);

% Predicciones del VAR
for i = 1:N
    yhat(:,i) = results(i).yhat;
end
P = eye(N); 
vnames =  ['IBEX           ',
           'SABADELL       ',    
           'BBVA           ',
           'SANTANDER      ',
           'INDICE_BANCARIO',
           'BANKIA         '];
names={'IBEX' 'SABADELL' 'BBVA' 'SANTANDER' 'INDICE BANCARIO' 'BANKIA'};

plt_var(results,vnames); 
prt_var(results,vnames);

% Ninguna es significativa salvo:
%
% - SAN: IBEX y BKA
% - IB : IBEX y BKA

%% Análisis de los residuos
figure;
subplot(4,3,1);
autocorr(err(:,1),15);
title('ACF: IBEX' );
subplot(4,3,2);
autocorr(err(:,2),15);
title('ACF: Sabadell' );
subplot(4,3,3);
autocorr(err(:,3),15);
title('ACF: BBVA');
subplot(4,3,4);
autocorr(err(:,4),15);
title('ACF: SANTANDER');
subplot(4,3,5);
autocorr(err(:,5),15);
title('ACF: ÍNDICE BANCARIO');
subplot(4,3,6);
autocorr(err(:,6),15);
title('ACF: Bankia');
subplot(4,3,7);
parcorr(err(:,1),15);
title('PACF: IBEX');
subplot(4,3,8);
parcorr(err(:,2),15);
title('PACF: SABADELL');
subplot(4,3,9);
parcorr(err(:,3),15);
title('PACF: BBVA');
subplot(4,3,10);
parcorr(err(:,4),15);
title('PACF: SANTANDER');
subplot(4,3,11);
parcorr(err(:,5),15);
title('PACF: ÍNDICE BANCARIO');
subplot(4,3,12);
parcorr(err(:,6),15);
title('PACF: BANKIA');

%% Analizamos conectividad mediante una indentificación generalizada.

[MFRI,MMA]=fri(results, MPHI, NAR, N, M, P, vnames);
DeltaC = zeros(N,N);
for i = 1:N
    for j=1:N
        num=0;
        den=0;
        for h=1:H
            ei=zeros(N,1); ei(i)=1;
            ej=zeros(N,1); ej(j)=1;
            num = num+(ei' * MMA(1+N*(h-1):N*h,:)*Omega*ej)^2;
            den = den+ei'*MMA(1+N*(h-1):N*h,:)*Omega*MMA(1+N*(h-1):N*h,:)'*ei;
        end
        DeltaC(i,j)=(1/Omega(j,j))*num/den;
    end
end

Deltat = zeros(N,N);
for i=1:N
    for j=1:N
        S=sum(DeltaC');
        Deltat(i,j)=DeltaC(i,j)/S(i);
    end
end

disp('_____________________________________________________');
disp('Horizonte de previsión en la descomposición de la varianza');
H
disp('_____________________________________________________');
Deltat=Deltat.*100;

disp('_____________________________________________________');
disp('Net pairwise directional connectedness measures');
disp('Cada elemento (i,j) mide el efecto neto: en qué porcentaje i explica a j menos en qué porcentaje j explica a i');
Deltatnet=zeros(N,N);
for i=1:N
    for j=1:N
        Deltatnet(i,j)=Deltat(j,i)-Deltat(i,j);
    end
end
Deltatnet
disp('_____________________________________________________');
disp('Total directional Connectedness from others to i');
disp('Mide qué porcentaje de la varianza del error de previsión de i es explicado por el resto de las variables del VAR');
sum(Deltat,2)-diag(Deltat)
disp('_____________________________________________________');
disp('Total directional Connectedness from j to others');
disp('Mide qué porcentaje de la varianza del error de previsión de todas las variables (excepto j) es explicado por j');
sum(Deltat,1)-diag(Deltat)'
disp('_____________________________________________________');
disp('Net total directional connectedness');
disp('Mide el porcentaje que explica i de los otros - el porcentaje de i explicado por otros');
(sum(Deltat,1)-diag(Deltat)')'-(sum(Deltat,2)-diag(Deltat))
disp('_____________________________________________________');
SS=0;
for i=1:N
    for j=1:N
        if i~=j
            SS=SS+Deltat(i,j);
        end
    end
end
disp('Total connectedness');
SS=SS/N
disp('_____________________________________________________');

%% HACEMOS GRÁFICO
ss=[];
tt=[];
weights=[];
for i=1:N
    for j=1:N
        if 100*Deltatnet(i,j)./max(Deltatnet)>20
            ss=[ss,i];
            tt=[tt,j];
            weights=[weights,Deltatnet(i,j)];
        end
    end
end
weights1 = round(weights/max(weights)*100);LWidth=4/100*weights1;
G = digraph(ss,tt,weights1,names);
%VNodeCData=[0 0.5 1 1.5 2 2.5 3 3.5 4 4.5];
VNodeCData=[0 0.5 1 1.5 2 2.5];
plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',LWidth,...
    'ArrowSize',10,'markersize',20,'NodeCData',VNodeCData,'EdgeColor','k')
axis equal off