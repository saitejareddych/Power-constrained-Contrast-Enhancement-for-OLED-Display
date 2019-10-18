clc;
close all
clear all
[J P]=uigetfile('*.*','select the source file');
I=imread(strcat(P,J));
I=imresize(I,[256 256]);
save I I
I=im2double(I);
Pow=.1;                                       % 10 percent
L=I.*sqrt(1-Pow);
%==============================
I=L;
P=sum(sum(I(:,:,1)))+sum(sum(I(:,:,2)))+sum(sum(I(:,:,3)));
% %===============================================
% Proposed Approach ======================
Ycr=rgb2ycbcr(I);
gamma=.5;
Tau=127;
D=256;
YP=sum(sum(Ycr(:,:,1).^gamma));
WL=((Tau/(D-1))*log(D))/(log(Tau+1));
WH=((1-(Tau/(D-1)))*log(D))/(log(D-Tau));
alp=.1;
for N=1
for ii=1:size(Ycr,1)
    for jj=1:size(Ycr,2)
        if Ycr(ii,jj,N)<=Tau
            nm=Tau*log(D)*log(alp*Ycr(ii,jj,N)+1);
            dm=(D-1)*log(alp*Tau+1);
            plog(ii,jj,N)=nm/dm;     
        else          
            mlog=-WH*log(D-Ycr(ii,jj,N))+log(D);
            plog(ii,jj,N)=mlog;            
        end
    end
end
end

[M N Z]=size(I);
[Y,X]=meshgrid(1:N,1:M);
 c=[50 170 220];
 n=length(c);
 RF=zeros(size(plog));
for ii=1:n 
 Fnok = exp(-((X.^2)+(Y.^2))./(c(n).^2));
 K = 1/(sum(sum(Fnok)));
 F = K.*Fnok; 
 IR = double(plog);
 FF  = fftshift(fft2(F));
 IFR = fftshift(fft2(IR));
 IFR=FF.*IFR;
 IFR=real(ifft2(ifftshift(IFR)));
 RR = double(IR)-(IFR);
 nrn=abs(RR)./max(max(RR));
 g=(1./(nrn+0.001)).^(1-c(n)/(max(c)+0.001));
 RF=RF+mat2gray(RR).*g;
end
yr=sum(sum(RF)).^gamma;
p=1-yr/YP;
fx=(plog-min(min(RF))) /(max(max(RF))-min(min(RF)));
del=(P.*Pow-p)/(M*N*Z); % for 10 %
T1=0.005;T2=0.1;
e=0.01;pp=.1;
[Rout nwtau]=newtt(Tau,T1,T2,del,e,Ycr(:,:,1),YP,P,pp);
R(:,:,1)=mat2gray(Rout);
R(:,:,2)=(Ycr(:,:,2));
R(:,:,3)=(Ycr(:,:,3));
RR=ycbcr2rgb(R);
g=1e-4:.5e-3:0.01;
figure,plot(g,log10(1./g));xlabel('---e');ylabel('Iterations')
ylim([1 5]);
figure,subplot(221);imshow(I);title('Original Image');
subplot(222);imshow(RR);title('proposed Enhanced Image');
disp('For Proposed ');
E1=sum(sum(20*log(max(RR)./(min(RR)+eps))))/(256*256)*10
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(RR(:,:,1)), hy, 'replicate');
Ix = imfilter(double(RR(:,:,1)), hx, 'replicate');
g1 = sqrt(Ix.^2 + Iy.^2);
G1=mean2((g1))*(Pow*100)


%==================================================
%==================================================
%clear 
load I
I=double(I);
P=.1;                                       % 10 percent
L=I.*sqrt(1-Pow);
subplot(223),imshow(uint8(L));title('Linear Image');
LL=mat2gray(L);
disp('For Linear ');
E2=sum(sum(20*log(max(LL)./(min(LL)+eps))))/(256*256)*10
I=L;
Iy = imfilter(double(LL(:,:,1)), hy, 'replicate');
Ix = imfilter(double(LL(:,:,1)), hx, 'replicate');
g2 = sqrt(Ix.^2 + Iy.^2);
G2=mean2((g2))*(Pow*100)
%=====================================
[Y,U,V] = rgb2yuv(I(:,:,1),I(:,:,2),I(:,:,3));
Yn = double(Y);
% extract an input histogram vector h
h = zeros(256,1);
for j=1:size(Yn,1)
    for i=1:size(Yn,2)
        temp = Yn(j,i) + 1;
        h(temp,1) = h(temp,1) + 1;
    end
end
clear temp
% modified histogram by logarithmic function
m = LHM(h,.5);
% Convex optimization
D = inv(tril(ones(256,256)));
[y, ~, ~] = PCCE(m, h, .1); x = D\y;
% write output image
for j=1:size(Yn,1)
    for i=1:size(Yn,2)
        out_Y(j,i) = round(x(Yn(j,i)+1));
    end
end
out_RGB = yuv2rgb(uint8(out_Y),U,V);
subplot(224),imshow(out_RGB);title('PCCE output');
out_RGB=double(out_RGB);
disp('For PCCE ');
E3=sum(sum(20*log(max(out_RGB)./(min(out_RGB)+eps))))/(256*256)*10
clear Ix Iy
Iy = imfilter(mat2gray(out_RGB(:,:,1)), hy, 'replicate');
Ix = imfilter(mat2gray(out_RGB(:,:,1)), hx, 'replicate');
g3 = sqrt(Ix.^2 + Iy.^2);
G3=mean2((g3))*(Pow*100)
