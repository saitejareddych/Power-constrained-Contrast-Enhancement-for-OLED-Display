clc;
clear all;
close all;
Tau=128;
D=256;
pix=[5:5:250];
alp=[1 .3 .1 0.03 .01];
cl={'--','r--','k--','m--','g--'};
WL=((Tau/(D-1))*log(D))/(log(Tau+1));
WH=((1-(Tau/(D-1)))*log(D))/(log(D-Tau));
for aa=1:length(alp)
    for bb=1:length(pix)
        if pix(bb)<=Tau
            nm=Tau*log(D)*log(alp(aa)*pix(bb)+1);
            dm=(D-1)*log(alp(aa)*Tau+1);
            plog(bb)=nm/dm;     
        else          
            mlog(bb)=-WH*log(D-pix(bb))+log(D);
            plog(bb)=mlog(bb);            
        end
    end
    plot(pix,(plog),cl{aa},'Linewidth',2);hold on;    
end
grid on;
xlabel('Pixel Values');
ylabel('--Plog');
legend('alp=1','alp=0.3','alp=0.1','alp=0.03','alp=0.01',4)


