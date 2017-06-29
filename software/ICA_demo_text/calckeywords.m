function str=calckeywords(T,A,terms,keywordstreshold,sortorder)
% Calc keywords

if nargin<5,
  sortorder=1:size(A,2);
end

dim=size(A,2);
try
  W=T(:,1:dim)*A;    
catch
  W=T*A;    
end
Wn=W;%Wn=W*diag(sign(mean(W))); % vende de negative
%keyindex=sign(mean(W))
Wn=Wn*inv(diag(max(Wn))); % normere til 1
str='';
is=0;
disp(sprintf('Keyword treshold %0.2f\n',keywordstreshold))
for i=sortorder
  is=is+1;
  %subplot(dim,1,is);
  %plot(Wn(:,i));
  indexkeywords=find(Wn(:,i)>keywordstreshold);
  [dummy,sortindex]=sort(Wn(indexkeywords,i));
  w=fliplr(Wn(indexkeywords(sortindex),i)');
  k=fliplr(terms(indexkeywords(sortindex)));
  keyword='Keywords: ';
  weight ='Weights : ';
  for j=1:length(w),
    keyword=[keyword,' ',char(k(j))];
    weight= [weight,' ',sprintf('%0.2f',w(j))];
  end
  disp(sprintf('Group %i',i))
  disp(keyword)
  disp(sprintf('%s\n',weight))
  str=[str sprintf('<B>Group %i</B><BR>',i) keyword sprintf('<BR>%s<P ALIGN=left><IMG SRC="timecomp%i.jpg"></P><BR>\n',weight,i)];
end
