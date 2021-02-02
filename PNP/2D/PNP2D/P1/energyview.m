clear
clc

oe = textread('energy.txt');
t = oe(3:end,1);
e = oe(3:end,2);

std=0*5*40*40*ones(length(e),1);
e=e-std;

plot(t,e)
title('P^2, N=64 \times 64, \Delta t = 0.001, T=60');
xlabel('t')
ylabel('Energy')
%axis([0 20 -5 60])
