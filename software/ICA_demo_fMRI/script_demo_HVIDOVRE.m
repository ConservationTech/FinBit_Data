function script_demo_HVIDOVRE

% Default settings
NrCmp=5;
ptreshold=0.025; % Image color treshold for % of pixels

% Load data
load dataHVIDOVRE.mat
X=X-repmat(mean(X')',1,size(X,2));
X=X/std(X(:));

imTranspose=0;
showMethode=1;

infofigure('HVIDOVRE fMRI Data','infoHVIDOVRE.txt');

CommandMenu=1;
while CommandMenu~=0
   CommandMenu=menu(sprintf('ICA [Sources:%d]',NrCmp)...
      ,'Sources  ->                  '...
      ,'Sources <-                   '...
      ,'Run: PCA                     '...
      ,'Run: MS-ICA                  '...
      ,'Run: ML-ICA                  '...
      ,'Run: MF-ICA                  '...
      ,'Run: MF-ICA (positiv sources)'...
      ,'[ Back ]                     ');
  
   switch CommandMenu,
      case 1, NrCmp=NrCmp+1;
      case 2, NrCmp=NrCmp-1;
      case 3, run_ica(0,[],NrCmp,X,P,dim,imTranspose,showMethode,ptreshold);
      case 4, run_ica(2,[],NrCmp,X,P,dim,imTranspose,showMethode,ptreshold);
      case 5, run_ica(1,[],NrCmp,X,P,dim,imTranspose,showMethode,ptreshold);
      case 6, run_ica(3,[],NrCmp,X,P,dim,imTranspose,showMethode,ptreshold);
      case 7, run_ica(4,[],NrCmp,X,P,dim,imTranspose,showMethode,ptreshold);
      otherwise, CommandMenu=0;
   end

end
close all