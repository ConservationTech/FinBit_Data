function [cM,cMrel,nrInClasses]=confusionMatrix(Targets,Estimats)
% confusionMatrix.m     : Calculate the confusion matrix where cols are the target 
%                         classes and rows are the estimats.
%
%                         [cM,cMrel,nrInClasses]=confusionMatrix(Targets,Estimats)
%
%                         In:
%                         Targets  : Vector of target classes
%                         Estimats : Vector of estimatet classes
%
%                         Out:
%                         cM : confusion matrix.
%                         cMrel : confusion matrix that sum to 100 procent in the cols.
%                         nrInClasses : Vector with number of results for each class
%
% (29.4.99 TK)

nrClasses=max(Targets);
nrEstimats=max(Estimats);

%cM=zeros(nrClasses,nrClasses);
for tc=1:nrClasses,
  Tindex=find(Targets==tc);
  for ec=1:nrEstimats,
    Eindex=find(Estimats(Tindex)==ec);
    cM(ec,tc)=length(Eindex);
  end
end

% Cols sum to 100 procent
cMrel=cM*diag(1./(sum(cM)))*100;

% sort rows
%p = zeros(1,nrClasses);
%c=cM;
%for i=1:nrClasses
%  [coeff,p(i)] = max(c(i,:));
%  c(:,p(i))=zeros(nrClasses,1)-1;
%end

%cMrel = cMrel(:,p);
%cM = cM(:,p);

nrInClasses=sum(cM);
