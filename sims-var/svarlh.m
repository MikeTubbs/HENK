function llh=svarlh(x,S,idmat,T)
% function llh=svarlh(x,S,idmat,T)
% x:  vector of free parameters
% S:  estimated cross product matrix (not moments) of rf residuals
%--------------------
nv=size(idmat,1);
%---------------------
aloc=find(idmat);       
A=zeros(4,4);
A(aloc)=x;
[u,d,v]=svd(A);
llh=sum(log(diag(d)))*T-.5*sum(sum((v*d.^2*v').*S));
llh=-llh; % reverse sign for use in minimization routine
