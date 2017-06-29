function infofigure(title,filename);

% Display intro text
h=figure(1); clf
set(h,'Name',title);
set(h,'NumberTitle','off');
colorf=get(h,'Color');
posf=get(h,'Position');

% Title
pos = [10 posf(4)-40 posf(3)-20 30];
h = uicontrol('Style','Text','Position',pos);
set(h,'String',{title},'FontSize',14)

% Info text
pos = [10 10 posf(3)-20 posf(4)-60];
h = uicontrol('Style','Text','Position',pos,'BackgroundColor',colorf);

i=0;
fid = fopen(filename);
while 1
    i=i+1;
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    string{i}=tline;
end
fclose(fid);

set(h,'String',string,'HorizontalAlignment','Left','FontSize',12)