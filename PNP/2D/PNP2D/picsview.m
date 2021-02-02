clear
close all
clc

Fig=0;
target=2;

FileNamestop='';
xFileNamestop='';

while 1 %& (Fig==target)
    FileName=strcat('ue', int2str(Fig), '.txt');
    xFileName=strcat('x.txt');

    fp=fopen(FileName);
    if (fp<0)
        break;
    end
    fclose(fp);
       
    Fig=Fig+1;
    Ox = textread(xFileName);
    Ou = textread(FileName);
    M=Ox(1,1);
    t=Ou(1);
    x=Ox(2:end,1);
    y=Ox(2:end,2);
    u=Ou(2:end);
    
for i=1:M
    for j=i:M
       if (x(j) < x(i))
           xtemp = x(i);
	       x(i) = x(j);
	       x(j) = xtemp;
	       
           ytemp = y(i);
	       y(i) = y(j);
	       y(j) = ytemp;
	       
	       utemp = u(i);
	       u(i) = u(j);
      	   u(j) = utemp;
        end
    end
end

i=1;
j=2;
while i<=M
     while j<=M
         if ( (x(j) ~= x(i) & x(j-1) == x(i)) | j==M )
           if (j==M)
               j=j+1;
           end
           [i j-1 x(i) x(j-1)];
           for ii=i:j-1
               for jj=ii:j-1
                 if (y(jj) < y(ii))
                   ytemp = y(ii);
                   y(ii) = y(jj);
                   y(jj) = ytemp;
                   
                   utemp = u(ii);
	               u(ii) = u(jj);
      	           u(jj) = utemp;
                 end
               end
           end
           i=j;
         end
         j=j+1;
     end
end

% The following is for plot only
i=1;
j=2;
ix = 1;
while i<=M
     while j<=M
         if ( (x(j) ~= x(i) & x(j-1) == x(i)) | j==M )
           if (j==M)
               j=j+1;
           end
           ii=i;
           jj=ii+1;
           while ii<=j-1
               while jj<=j-1
                   if (y(jj-1) == y(ii) & ii==j-1)
                       xx(ix) = x(i);
                       yy(ix) = y(ii);
                       uu(ix) = u(ii);
                       ix=ix+1;
                       ii=jj+1;
                       disp('I am here')
                   end
                   if (y(jj) ~= y(ii) & y(jj-1) == y(ii))
                       xx(ix) = x(i);
                       yy(ix) = y(ii);
                       u_sum=0.0;
                       for k=ii:jj-1
                           u_sum=u_sum+u(k);
                       end
                       uu(ix) = u_sum/(jj-ii);
                      % jj-ii
                       ix=ix+1;
                       ii=jj;
                   end
                   if (y(jj) == y(ii) & jj==j-1)
                       xx(ix) = x(i);
                       yy(ix) = y(ii);
                       u_sum=0.0;
                       for k=ii:jj
                           u_sum=u_sum+u(k);
                       end
                       uu(ix) = u_sum/(jj-ii+1);
                       ix=ix+1;
                       ii=jj+1;
                   end
                   jj=jj+1;
               end
           end
           i=j;
         end
         j=j+1;
     end
end

xx=xx';
yy=yy';
uu=uu';

MM=sqrt(length(uu));

uuu=zeros(MM,MM);
[xxx,yyy]=meshgrid(yy(1:MM));
for i=1:MM
    for j=1:MM
        uuu(j,i)=uu( MM*(i-1)+j );
    end
end
% uuu

isfig = 1;

if (isfig~=0)
    
figure(Fig)
surf(xxx,yyy,uuu)

figure(Fig+20)

colormap(jet(128))

 p=pcolor(xxx,yyy,uuu)
 colorbar
 title('Pattern formation');
 set(p,'EdgeColor', 'flat');

%title(['Pattern formation wtih P^3, N=20, at t=', num2str(t)]);
%xlabel('x')

%ylabel('u_h')
%axis([0 8*pi -4 4])
%axis([0 8*pi -1.0 1.0])
hold on

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


