function [S,A]=unitvariance(S,A)
% Variance = 1 on outputs - corrected in A matrix

A=A.*((diag(std(S'))*ones(size(A,1)))');
S=S.*(inv(diag(std(S')))*ones(size(S,1),size(S,2))); 
