function [Rout tau]=newtt(tau,T1,T2,del,e,Ycr,YP,P,pp)
ta=tau;lamda=0;
it=5;
for tt=1:it
%if del>-T1 && del<T1      
if abs(del)<e
    Rout=Ycr;
else
   
    if -T2<del && del<-T1
        ta=tau+(D-tau)/4;
        tau=ta;
    end
    
    if -T2<del && del<-T1
        ta=tau+(D-tau)/4;
        tau=ta;
    end
    
    if T1<del && del<T2
        ta=tau-tau/4;
        tau=ta;
    end
    
    if del<-T2
    ta=tau+(D-tau)/2;
    tau=ta;
    end
    
    if del>T2
    ta=tau-tau/2;
    tau=ta;
    end
    
    [dd1 Ycr]=algo(tau,Ycr,YP,P,pp);    
    del=dd1;
    Rout=Ycr+lamda*Ycr;
    lamda=lamda-del;
end

end
%==================================================
function [del plog]=algo(Tau,Ycr,YP,P,pp)
D=256;gamma=0.5;
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
M =256;
N =256;
Z=3;
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
del=(P.*pp-p)/(M*N*Z); % for 10 %    
       

    
    