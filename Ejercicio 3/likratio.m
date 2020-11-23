function likratio(y,p1,p0,prob,simsc)
% PURPOSE: performs likelihood ratio test for var model
%---------------------------------------------------
% USAGE:  likratio(y,p1,p0,simsc) 
% where:    y     = an (nobs x neqs) matrix of y-vectors
%           p1    = retardos del VAR sin restringir
%           p0    = retardos del VAR restringido (p0<p1)
%           prob  = probabilidad de la significación del contraste
%           simsc = flag for Sim's dof correction factor 
%                    0 = no, 1 = use correction
%                    (default = 0)
%---------------------------------------------------
% RETURNS: nothing, prints results to the MATLAB command window
%---------------------------------------------------
% SEE ALSO: var, varf, prt_var 
%---------------------------------------------------

if nargin > 5
error('wrong # of arguments to lrratio');
elseif nargin == 4
simsc = 0;
end;

if p1 < p0
 error('p1 < p0 in likratio');
end;

if p0  < 1
    p0 = 1;
end;

[nobs n] = size(y);

resid1=var_resid(y,p1);
resid0=var_resid(y,p0);
omega1=cov(resid1);
omega0=cov(resid0);
T=nobs-p1;
if simsc==1
    c=n*p1-1;
else
    c=0;
end
LR=(T-c)*(log(det(omega0))-log(det(omega1)))
vc=chi2inv(1-prob,(n^2)*(p1-p0))
pvalue=1-chi2cdf(LR,(n^2)*(p1-p0))

