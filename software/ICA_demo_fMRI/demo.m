function demo
% demo.m            Demonstration of fMRI analysis with Indepandent component analysis (ICA)

% by Thomas Kolenda 2002

format compact

CommandMenu=1;
while CommandMenu~=0
   CommandMenu=menu('ICA - fMRI demos'...
      ,'MAPAWAMO demo'...
      ,'HVIDOVRE demo'...
      ,'[ End ]                      ');
  
   switch CommandMenu,
      case 1, script_demo_MAPAWAMO;
      case 2, script_demo_HVIDOVRE;
      otherwise, CommandMenu=0;
   end

end