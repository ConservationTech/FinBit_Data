function run_ica(algorithm,slice,nrcmp,XN,P,dim,imTranspose,showMethode,ptreshold)

figure(1)
clf
title('Working...')
axis off
drawnow

try
    slice=slice;
catch
    slice=1;
end    

% rescale data
XN=XN/max(max(XN));

P=(P-.5)*1.8;
M=nrcmp;

svdT=0;
icaT=0;

if algorithm==0
    tic
    [U,D,V]=svd(XN',0);
    svdT=toc;
    S=U(:,1:M)';
    T=V(:,1:M)*D(1:M,1:M);
    S=(U(:,1:M)*D(1:M,1:M))';
    T=V(:,1:M);
    str=['PCA'];    
end

if algorithm==1
    cd icaML
    tic
    [U,D,V]=svd(XN',0);
    svdT=toc;
    tic
    [S,A]=icaML((U(:,1:M)*D(1:M,1:M))');
    icaT=toc;        
    T=V(:,1:M)*A;
    [S,T]=ica_sortnflip(S,T);
    str=['ICA ML'];
    cd ..
end

if algorithm==11
    cd icaML    
    tic
    [U,D,V]=svd(XN',0);
    svdT=toc;    
    tic
    [S,A]=icaML(D(1:M,1:M)*V(:,1:M)');
    icaT=toc;            
    T=S';
    S=(U(:,1:M)*A)';
    [S,T]=ica_sortnflip(S,T);
    str=['ICA ML time'];
    cd ..
end

if algorithm==12
    cd icaMS    
    tic
    [U,D,V]=svd(XN',0);
    svdT=toc;    
    tic
    [S,A]=icaMS((U(:,1:M)*D(1:M,1:M))');
    icaT=toc;            
    T=V(:,1:M)*A;
    [S,T]=ica_sortnflip(S,T);
    str=['ICA MS'];
    cd ..
end

if algorithm==2
    cd icaMS    
    tic
    [U,D,V]=svd(XN',0);
    svdT=toc;    
    tic
    [S,A]=icaMS(D(1:M,1:M)*V(:,1:M)');
    icaT=toc;            
    T=S';
    S=(U(:,1:M)*A)';
    [S,T]=ica_sortnflip(S,T);
    str=['ICA MS time'];
    cd ..
end

if algorithm==3 | algorithm==4
    cd icaMF    
    rand('seed',0); randn('seed',0);
    if algorithm==3 
        prior.S='heavy_tail'; 
        str=['ICA MF'];
    else
        prior.S='exponential'; 
        str=['ICA MF pos sources'];
    end
    par.sources=M;par.S_max_ite=20;
    tic
    [S,A]=ica_adatap(XN,prior,par,1);
    S=real(S); A=real(A);
    icaT=toc;            
    [S,A]=ica_sortnflip(S,A);    
    T=A;
    cd ..
end

figure(1)
clf
col=4;

Xbackground=reshape(mean(XN),dim(1),dim(2));
Xbackground=((Xbackground-min(Xbackground(:)))/max(Xbackground(:)-min(Xbackground(:)))*63)+1;

if imTranspose==3,
    Xbackground=flipud(Xbackground');
end

for i=1:M
    bg=Xbackground;
    % IC image
    subplot(M,col,i*col-col+1)
    x=reshape(S(i,:),dim(1),dim(2));
    x=(x)/max( abs(max(x(:))) , abs(min(x(:))) )*32+32;
    
    if imTranspose==3,
        x=flipud(x');
    end
    
    xs=sort(x(:));
    pu=round(xs(round(length(x(:))*(1-ptreshold))));
    pl=round(xs(round(length(x(:))*(ptreshold))));

    if showMethode==1,
        bg(find(x>pu))=64;
        bg(find(x<pl))=1;
        
        co=gray;
        co(1,:) =[0,0,0.5625];
        co(64,:)=[0.5625,0,0];
        
        figure(1);
        image(bg);
        colormap(co);
        
    elseif showMethode==2
        g=gray; j=jet;
        if algorithm==4,
            co=[g(1:pu-1,:);j(pu:end,:)];
        else
            co=[j(1:pl,:);g(pl+1:pu-1,:);j(pu:end,:)];
        end
        
        figure(1);
        image(x);
        colormap(co);        
    end

    set(gca,'YTicklabel','');
    set(gca,'XTicklabel','');  
    axis square
    
    % IC time
    subplot(M,col,i*col-2:i*col)
    Ti=T(:,i)/max(abs(T(:,i)));
    plot(1:length(P),Ti,'-',1:length(P),P,':');
    axis([1 length(P) -1 1])
    set(gca,'YTick',[]);
    set(gca,'XTick',[]); 
    
    if i==1
        tit=sprintf('%s',str);
        title(tit)
    end;
end
drawnow

% Display result text
h=figure(1);
posf=get(h,'Position');

str='';
if svdT~=0
  str=[str sprintf('PCA time:%0.2f  ',svdT)];
end
if icaT~=0
  str=[str sprintf('ICA time:%0.2f  ',icaT)];
end
str=[str 'seconds'];

pos = [10 10 posf(3)-20 30];
h = uicontrol('Style','Text','Position',pos);
set(h,'String',{str},'FontSize',10)


% try
%    print('-depsc',['./',tit]);
%    save(tit)
% catch
% end


function [S,A]=ica_sortnflip(S,A);
% Sort and assign sign components
[M N] = size(S);

% Sort according to variance
Avar=diag(A'*A)/M;
Svar=diag(S*S')/N;
sig=Avar.*Svar;
[a,indx]=sort(sig);
S=S(indx(M:-1:1),:);
A=A(:,indx(M:-1:1));

% flip components w.r.t. having positive mean
Per=sign(mean(S'));
S=diag(Per)*S;
A=A*diag(Per);


