function plotComp(X,base,header,b)

DS = X';

if nargin<2,
   base=0;
   header='IC comps';
end

try
   if base==0,
      base=ones(1,size(DS,1));
   end
end

clf
if isstruct(base)==1,
   if min(size(cell2mat(base.qrels')))>1,
      [non,group]=max(cell2mat(base.qrels')');
      num_groups=min(size(cell2mat(base.qrels')))
   else
      group=cell2mat(base.qrels');
   end
   
else
   group=base;
end

maxNumGroups=max(group);
maxGroups=[];
for i=1:maxNumGroups,
   if sum(i==group)>0,
      maxGroups=[maxGroups i];
   end
end
maxGroups=length(maxGroups);

maxComp=size(DS,2);
if maxComp>5, 
   maxComp=5;
end

color='rgbkcmy';


if nargin<4,
   b=2;
end

a=[];
for i=1:5,
   if i~=b,
      a=[a i];
   end
end

if maxGroups<4,
   showGroup=maxGroups;
else
   showGroup=4;
end
maxComp;

l='';
for n=1:maxComp-1
   subplot(2,2,n)
   if maxGroups>length(color),
      plot(DS(:,a(n)),DS(:,b),'.b');
   else
      hold on
      for i=1:maxNumGroups,
         index=find(group==i);
         plot(DS(index,a(n)),DS(index,b),['.' color(i)]);
         if n==1,
            l=[l;mat2str(i)];
            title(header);
         end
      end
      hold off
   end
   xlabel(sprintf('COMP %i',a(n)))
   ylabel(sprintf('COMP %i',b))
end

if maxGroups<=length(color),
   subplot(2,2,1)
   %  legend(l)
end
