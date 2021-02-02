clear
clc

oe = textread('energy.txt');
t = oe(2:end,1);
e = oe(2:end,2);

plot(t,e)
title('Free energy with P^3, N=20');
xlabel('t')
ylabel('E')
%axis([0 20 -5 60])
