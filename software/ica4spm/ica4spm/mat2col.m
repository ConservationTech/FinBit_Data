function mat2col(im,back,frac)
[N,M]=size(im);
x=reshape(im,N*M,1);
bx=reshape(back,N*M,1);
[q,k]=sort(x);
bmin=min(bx);bmax=max(bx);
y=5+40*(bx -bmin)/(bmax-bmin);
indx_low=k(1:floor(N*M*frac(1)));
indx_high=k(ceil(N*M*frac(2)):(N*M));
y(indx_low)=ones(length(indx_low),1);
y(indx_high)=ones(length(indx_high),1)*64;
y=reshape(y,N,M);
image(y)
my_gray=gray;
my_gray(1,:)=[0 0 1];
my_gray(64,:)=[1 0 0];
colormap(my_gray)