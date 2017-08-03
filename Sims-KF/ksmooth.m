function [btT,StT]=ksmooth(btt,Stt,bnT,SnT,A,omega)
%[btT StT]=ksmooth(btt,Stt,bnT,SnT,A,omega)
% Smoothing recursion.  State evolution equation is
%    bn=A*bt+e,  Var(e)=omega
%    bt|t ~ N(btt,Stt) -- from Kalman Filter
%    bn|T ~ N(bnT,SnT) -- distribution of bn given full sample. From 
%                         KF if n=T, otherwise from this recursion
%    bt|T ~ N(btT,StT)
% btt and bnT are assumed to be loaded as row vectors and btT is also in
% the same format

% whos A Stt omega
AS=A*Stt;
G=AS*A'+omega;
%****************************
% SAGI=AS'/G;
% line below may slightly slow the routine, but makes it robust vs.
% cases where part or all of the state vector is known with certainty.
SAGI=AS'*pinv(G);
%*****************************
btT=(SAGI*(bnT'-A*btt'))'+btt;
StT=Stt-SAGI*AS+SAGI*SnT*SAGI';
