clear all; close all;
fs = 11025;
L = 2;
[xn, fs] = audioread('beethoven5.wav');
load('LPFlab1_2.mat')
k = 0:(length(LPF)-1);
HPF = (-1).^k.*LPF;
figure(1)
freqz(HPF);
figure(2)
freqz(LPF);
N=1000;
figure(3)
plotMagPhase(FTD((0:(length(xn)-1)),xn,N),N);

%% Option 1
%playsound(xn,fs);

%% Option 2
y2n = xn(1:L:end);
%playsound(y2n,fs/L);

%% Option 3
vn = conv(xn,LPF);

%playsound(vn,fs);
y3n = vn(1:L:end);
%playsound(y3n,fs/L);

%% Option 4
vn = conv(xn,HPF);
%playsound(vn,fs);
y4n = vn(1:L:end);
%playsound(y4n,fs/L);

%% Assignment 3: Signal reconstruction
w3n(1:L:length(y3n)*2) = y3n;
w3n = conv(w3n,LPF);

w4n(1:L:length(y4n)*2) = y4n;
w4n = conv(w4n,HPF);

yrec = 2*w3n + 2*w4n;
%playsound(yrec,fs);
figure(4);
plotMagPhase(FTD((0:(length(yrec)-1)),yrec,N),N);

figure(5);

plotMagPhase(FTD((0:(length(HPF)-1)),HPF+LPF,N),N);



function [  ] = playsound(x,fs)
    sound(x,fs);
    pause(length(x)*1/fs+.5);
end