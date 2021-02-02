clear
clc

Fig=0;
target=0;

% Reading the reference solution.

    Ox = textread('xxR.txt');
    Ou = textread('uuR.txt');
    MR=Ox(1,1);
    tR=Ou(1);
    xR=Ox(2:end,1);
    yR=Ox(2:end,2);
    uR=Ou(2:end);


while 1 % & (Fig==target)
    FileName=strcat('uu', int2str(Fig), '.txt');
    xFileName=strcat('xxR', '.txt');
    fp=fopen(FileName);
    if fp<0
        break;
    end
    fclose(fp);
    Fig=Fig+1;
    Ox = textread('xxR.txt');
    Ou = textread(FileName);
    M=Ox(1,1);
    t=Ou(1);
    x=Ox(2:end,1);
    y=Ox(2:end,2);
    u=Ou(2:end);
    
%    FileName=strcat('uR.txt');
%     fp=fopen(FileName);
%     if fp<0
%         break;
%     end
%     fclose(fp);
%    Fig=Fig+1;

    
 LinfErr=0.0;
 
 errcount = 1;
 for i=1:M
     for j=1:MR
         if(x(i)==xR(j) & y(i)==yR(j))
             if abs( u(i)-uR(j) )>=LinfErr
%             [i j x(i) y(i) xR(j) yR(j) abs( u(i)-uR(j) )]
             end
             LinfErr = max(  LinfErr , abs( u(i)-uR(j) ) );
             errcount = errcount +1;
         end
     end
 end

errcount;
LinfErr
    
    
    
    
    
    
    
    
isfig =0;

if (isfig~=0)
    
figure(Fig)
%surf(xxx,yyy,uuu)
plot3(x,y,u,'.')
title(['Pattern formation wtih P^3, N=20, at t=', num2str(t)]);
%xlabel('x')
%ylabel('u_h')
%axis([0 8*pi -4 4])
%axis([0 8*pi -1.0 1.0])
%hold on

% figure(6)
% if Fig~=1
% plot(x,u)
% end
% title('Pattern formation wtih P^3, N=20 for t > 0.0');
% xlabel('x')
% ylabel('u_h')
% axis([0 8*pi -1.0 1.0])
% hold on
end
 
 end


