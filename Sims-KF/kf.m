function [shatnew,signew,lh,yhat]=kf(y,H,shat,sig,G,M)
%function [shatnew,signew,lh,yhat]=kf(y,H,shat,sig,G,M)
% s is the state, and the plant equation is s(t)=Gs(t-1)+Me, where e is
% N(0,I).  The observation equation is y(t)=Hs(t).  The prior distribution for
% s is N(shat,sig).  To handle the standard Kalman Filter setup where the observation
% equation is y(t)=Hs(t)+Nu(t), expand s to become [s;v] (where v=Nu), expand H to [H I], replace
% G by [G 0;0 0], and replace M with [M 0;0 N].  The resulting model has no "error" in the
% observation equation but is equivalent to the original model with error in the state equation.
% The posterior on the new state is N(shatnew,signew) and lh is a two-dimensional vector containing
% the increments to the two component terms of the log likelihood function.  They are added 
% to form the log likelihood, but are used separately in constructing a concentrated or marginal
% likelihood. yhat is the forecast error of y based on the prior.
lh=zeros(1,2);
omega=G*sig*G'+M*M';
[uo do vo]=svd(omega);
[u d v]=svd(H*uo*sqrt(do));
first0=min(find(diag(d)<1e-12));
if isempty(first0),first0=min(size(H))+1;end
u=u(:,1:first0-1);
v=v(:,1:first0-1);
d=diag(d);d=diag(d(1:first0-1));
fac=vo*sqrt(do);
yhat=y-H*G*shat;
ferr=(v/d)*u'*yhat;
lh(1)=-.5*ferr'*ferr;
lh(2)=-sum(log(diag(d)));
shatnew=fac*ferr+G*shat;
signew=fac*(eye(size(v,1))-v*v')*fac';

