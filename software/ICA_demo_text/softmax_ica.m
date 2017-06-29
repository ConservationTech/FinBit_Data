function [S,A]=softmax_ica(S,A)

[M N] = size(S);

% vende de negative
[v,index]=max([abs(min(S'))' max(S')']');
S=diag(sign(index-1.5))*S;
A=A*diag(sign(index-1.5));

expS = exp(S);
S = expS./(ones(M,1)*sum(expS));

