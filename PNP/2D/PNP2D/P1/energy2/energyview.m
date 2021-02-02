close all
clear
clc

Fig=1;
target=1;
while 1 % & (Fig==target)
    FileName=strcat('energy', int2str(Fig), '.txt');
    fp=fopen(FileName);
    if fp<0
        break;
    end
    fclose(fp);
    
    oe = textread(FileName);
    N=41;
    edf=5*40*40*ones(N,1);
    if Fig==1
        t1 = oe(2:end,1);
        e1 = oe(2:end,2);
        step = (length(t1)-1)/(N-1);
        nt1 = zeros(N,1);
        ne1=nt1;
        for i=1:N
            nt1(i)=t1(step*(i-1)+1);
            ne1(i)=e1(step*(i-1)+1);
        end
        ne1=ne1-edf;
    end
    if Fig==2
        t2 = oe(2:end,1);
        e2 = oe(2:end,2);
        step = (length(t2)-1)/(N-1);
        nt2 = zeros(N,1);
        ne2=nt2;
        for i=1:N
            nt2(i)=t2(step*(i-1)+1);
            ne2(i)=e2(step*(i-1)+1);
        end
        ne2=ne2-edf;
    end    
    if Fig==3
        t3 = oe(2:end,1);
        e3 = oe(2:end,2);
        step = (length(t3)-1)/(N-1);
        nt3 = zeros(N,1);
        ne3=nt3;
        for i=1:N
            nt3(i)=t3(step*(i-1)+1);
            ne3(i)=e3(step*(i-1)+1);
        end
        ne3=ne3-edf;
    end    
    if Fig==4
        t4 = oe(2:end,1);
        e4 = oe(2:end,2);
        step = (length(t4)-1)/(N-1);
        nt4 = zeros(N,1);
        ne4=nt4;
        for i=1:N
            nt4(i)=t4(step*(i-1)+1);
            ne4(i)=e4(step*(i-1)+1);
        end
        ne4=ne4-edf;
    end    
    if Fig==5
        t5 = oe(2:end,1);
        e5 = oe(2:end,2);
        step = (length(t5)-1)/(N-1);
        nt5 = zeros(N,1);
        ne5=nt5;
        for i=1:N
            nt5(i)=t5(step*(i-1)+1);
            ne5(i)=e5(step*(i-1)+1);
        end
        ne5=ne5-edf;
    end
        
    Fig=Fig+1;
 
end    
plot(nt1,ne1,'-*',nt2,ne2,'-+',nt3,ne3,'-o',nt4,ne4,'-.',nt5,ne5,'-')
title('P^2, N=64 \times 64, T=10');
xlabel('t')
ylabel('Energy')
legend('\Delta t = 0.1','\Delta t = 0.1*2^{-1}','\Delta t = 0.1*2^{-2}','\Delta t = 0.1*2^{-3}','\Delta t = 0.1*2^{-4}' )
% axis([0 8*pi -1.2 1.2])
