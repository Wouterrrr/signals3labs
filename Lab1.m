%% Assignment 1: Sine Wave
f0 = 1.8e3; fs = 6e3;
n = 0:1000;
xn = sin(2*pi*f0*n./fs+pi/2);
load('lowpassfilter.mat');
fc = pi/3/2/pi*fs;
plotMagPhase(FTD((0:(length(lowPassFilter)-1)),lowPassFilter,N),N)