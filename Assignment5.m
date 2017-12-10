clear all; close all;

k = 3; M = 8;

N = 1e4;
n = -N:N;
an = sinc(k/M-n);
fig  = figure;
plotMagPhase(FTD(n,real(an),N),N);
saveas(fig,'an.png');
bn = (n-k==0);
fig = figure;

plotMagPhase(FTD(n,real(bn),N),N);
saveas(fig,'bn.png');