close all
clear
clc

x=-1:0.01:1;
y=ones(size(x));
plot(x,-1.*y,'k','LineWidth',1)
hold on
plot(-1.*y,x,'k','LineWidth',1)
plot(x,1.*y,'k','LineWidth',1)
plot(1.*y,x,'k','LineWidth',1)

x1 = [-1 0 1];
y1=ones(size(x1))

plot(x1,0.*y1,'*k','LineWidth',4)
plot(x1,-1.*y1,'*k','LineWidth',4)
plot(x1,1.*y1,'*k','LineWidth',4)

axis([-2,2,-1,3])
box off
axis off

set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
