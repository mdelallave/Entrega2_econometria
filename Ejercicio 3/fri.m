function [Mfri,MMA]=fri(results,MPHI,NAR,N,M,P,Vnames)
% results es el output de la función fvar.m de Le Sage
% MPHI son las matrices del VAR
% NAR es el número de retardos del VAR
% N es el número de variables del VAR
% M es el orden del VMA(inf). Por ejemplo 50.
% P es la matriz de factorización para calcular la fri de un modelo estuctural
% vnames es opcional: vector de nombres

% ----If no variable names supplied --------------------------
if nargin  == 3 
  Vnames = [];  
  for i=1:N
    Vnames = strvcat(Vnames,['Y' int2str(i)]);
  end
end
MMA=vma(MPHI,NAR,N,M);
MMAt=MMA*P';
friN=[];
for j=1:N
    for i=1:M
        frij(:,i)=MMAt((i-1)*N+1:i*N,j);
    end
    friN=[friN;frij];
end
x=[1:M]';
Mfri=[];
plotct=0;
for i=1:N;
    plotct = plotct + 1;
    plotdata = friN((i-1)*N+1:i*N,:)';
    Mfri=[Mfri plotdata];
    figure;
    plot(x,plotdata);
    title('funcion de respuesta a un impulso')
    ylabel(Vnames(plotct,:));
    xlabel(['Response of all variables to shock in equation ' num2str(plotct)]);    
    legend(Vnames);
    if plotct == N, break,  end
end
