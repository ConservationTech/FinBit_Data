function [ClassEstimats,ClassEstimatsRec]=classifyer_angel_ica(S,Frac);
% Classify the output of ICA using the largest value for a given document
% as class.

if nargin<2,
   Frac=0;
end

AllZeros=0;
LessDiff=0;
for i=1:size(S,2),
   [s,si]=sort(S(:,i));
   s=flipud(s);
   si=flipud(si);
   if s(1)==0,
      ClassEstimats(i,1)=size(S,1)+1;
      AllZeros=AllZeros+1;
   else
      ClassEstimats(i,1)=si(1);
      
      % reject class
      for j=2:length(s),
         %if s(j)>=s(1)-Frac,                               % NB! for softmax
         if abs(s(1)-s(j))/abs(s(1)+s(j))<Frac,
            ClassEstimats(i,1)=length(s)+1; 
            LessDiff=LessDiff+1;
            ClassEstimatsRec(i,j)=si(j);
         end
      end
   end
   
end

ClassEstimats=ClassEstimats';

disp(sprintf('Class - all zeros:%d    less frac diff:%d    = rejected:%d',AllZeros,LessDiff,AllZeros+LessDiff))
