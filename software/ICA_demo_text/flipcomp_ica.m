function [S,A]=flipcomp_ica(S,A)
% vende de negative IC componenter

[M N] = size(S);

[v,index]=max([abs(min(S'))' max(S')']');
S=diag(sign(index-1.5))*S;
A=A*diag(sign(index-1.5));


