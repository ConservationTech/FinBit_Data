function  script_demo_MAPAWAMO

% Default settings
Slice=2; 
SliceArray=[5,11];
NrCmp=5;
ptreshold=0.05; % Image color treshold for % of pixels

% Load data
load dataMAPAWAMO.mat
eval(sprintf('XN=X%d;',SliceArray(Slice)));
P=([zeros(1,20) ones(1,20) zeros(1,20) ones(1,20)]);
dim=[29,33];

imTranspose=3;
showMethode=2;

infofigure('MAPAWAMO fMRI data','infoMAPAWAMO.txt');

CommandMenu=1;
while CommandMenu~=0
   CommandMenu=menu(sprintf('ICA [Slice:%d - Sources:%d]',SliceArray(Slice),NrCmp)...
      ,'Slice  ->                    '...
      ,'Slice <-                     '...
      ,'Sources  ->                  '...
      ,'Sources <-                   '...
      ,'Run: PCA                     '...
      ,'Run: MS-ICA                  '...
      ,'Run: ML-ICA                  '...
      ,'Run: MF-ICA                  '...
      ,'Run: MF-ICA (positiv sources)'...
      ,'[ Back ]                     ');
  
   switch CommandMenu,
      case 1, if Slice==length(SliceArray),Slice=1;else Slice=Slice+1; end;
      case 2, if Slice==1,Slice=length(SliceArray);else Slice=Slice-1; end;
      case 3, NrCmp=NrCmp+1;
      case 4, NrCmp=NrCmp-1;
      case 5, run_ica(0,SliceArray(Slice),NrCmp,XN,P,dim,imTranspose,showMethode,ptreshold);
      case 6, run_ica(2,SliceArray(Slice),NrCmp,XN,P,dim,imTranspose,showMethode,ptreshold);
      case 7, run_ica(1,SliceArray(Slice),NrCmp,XN,P,dim,imTranspose,showMethode,ptreshold);
      case 8, run_ica(3,SliceArray(Slice),NrCmp,XN,P,dim,imTranspose,showMethode,ptreshold);
      case 9, run_ica(4,SliceArray(Slice),NrCmp,XN,P,dim,imTranspose,showMethode,ptreshold);
      otherwise, CommandMenu=0;
   end
   if CommandMenu == 1 | CommandMenu==2,
       eval(sprintf('XN=X%d;',SliceArray(Slice)));
   end;

end
close all