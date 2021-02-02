close all
clear
clc

x=0:0.01:3;
y=ones(size(x));
plot(x,y,'k','LineWidth',2)
hold on
plot(y,x,'k','LineWidth',2)
plot(x,0.*y,'k','LineWidth',2)
plot(0.*y,x,'k','LineWidth',2)
plot(x,2.*y,'k','LineWidth',2)
plot(2.*y,x,'k','LineWidth',2)
plot(x,3.*y,'k','LineWidth',2)
plot(3.*y,x,'k','LineWidth',2)
annotation('arrow',[0.39 0.65],[0.382 0.382],'LineWidth',2.5,'LineStyle','-','color',[0 0 1],'HeadStyle','plain');
annotation('arrow',[0.39 0.65],[0.653 0.653],'LineWidth',2.5,'LineStyle','-','color',[0 0 1],'HeadStyle','plain');
annotation('arrow',[0.389 0.389],[0.382 0.653],'LineWidth',2.5,'LineStyle','-','color',[0 0 1],'HeadStyle','plain');
annotation('arrow',[0.648 0.648],[0.382 0.653],'LineWidth',2.5,'LineStyle','-','color',[0 0 1],'HeadStyle','plain');

set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])

gtext('K_i');
gtext('K_{i_1}');
gtext('K_{i_2}');
gtext('K_{i_3}');
gtext('K_{i_4}');

gtext('e_{i_1}');
gtext('e_{i_2}');
gtext('e_{i_3}');
gtext('e_{i_4}');
