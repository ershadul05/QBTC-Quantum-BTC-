clc
clear all 
close all

%[Inputimage,pathname1]=uigetfile('*.*','Plz Select the cover image to be watermarked.');
x=imread('P_gray.jpg');
[row1 col1]=size(x);
x2=x(1:row1,1:col1);
[row col]=size(x2);
% x2=rgb2gray(x2);
x=x2;
x=double(x);
%x1=x;
y1=size(x); % matrix size
%m1=y1(1); %return the number of rows
%n1=y1(2); 
 %s=0;
 BR=0;
 blocksize=16;
 Z=zeros(row1, col1);
for i=1:blocksize:row
for j=1:blocksize:col
    y=x(i:i+(blocksize-1),j:j+(blocksize-1));
    m= mean(mean(y));
    sig=std2(y);
    b=y>m;
    %t([i:i+(blocksize-1)],[j:j+(blocksize-1)])=b([1:blocksize],[1:blocksize]);
    k=sum(sum(b));
    if (k~=blocksize^2)&(k~=0)
    ml=m-sig*sqrt(k/((blocksize^2)-k));
    mu=m+sig*sqrt(((blocksize^2)-k)/k);
    mlr=abs(round(ml));
    mur=abs(round(mu));
    mlrd2b=dec2bin(mlr,8);
    
    mlrd2b=uint16(mlrd2b)-48;
    ml_n_one=numel(mlrd2b);
    murd2b=dec2bin(mur,8);
    murd2b=uint16(murd2b)-48;
    mu_n_one=numel(murd2b);
    
    [u v w]=find(b);
    sbit=nnz(w);
    rowu=dec2bin(u,log2(blocksize));
    urowu=uint16(rowu);
    aurowu=urowu-48;
    erowu=numel(aurowu);
    %erowuall=numel(erowu);
    
    rowv=dec2bin(v,log2(blocksize));
    urowv=uint16(rowv);
    aurowv=urowv-48;
    erowv=numel(aurowv);
    %ecolvall=numel(aurowv);
    
    statebit_one=erowu+erowv;
    
    numone_b=nnz(w);
    
    threshold_bits=ml_n_one+mu_n_one;
    
    BR=BR+(numone_b+(statebit_one)+threshold_bits+sbit)/(1024*1024);
    
    x(i:i+(blocksize-1),j:j+(blocksize-1))=b*mu+(1-b)*ml;

    Z(i:i+blocksize-1,j:j+blocksize-1)=b;
end
end
end
imshow(x,[0 255])
%imshow(Z);
numofblockr=row1/blocksize;
numofblockc=col1/blocksize;

rposioflastblock=row1-blocksize;
cposioflastblock=col1-blocksize;

numofcolbit=numel(dec2bin(rposioflastblock))-numel(dec2bin(blocksize));
numofrowbit=numel(dec2bin(cposioflastblock))-numel(dec2bin(blocksize));
tbitr=(numofblockr*numofblockc*(numofcolbit+numofrowbit))/(1024*1024);

BR=BR+tbitr;
%figure(1);
%subplot(131),imshow(x2),title('ORIGINAL');
%subplot(132),imshow(t);title('ENCODED');
%subplot(133),imshow(uint8(x));title('DECODED');
BPP=(BR*1024*1024)/(row1*col1)
PSNR1=CalculatePSNR(x2,x)