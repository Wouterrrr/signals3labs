clear all; close all;

%% Assignment 1: Sine Wave
f0 = 1.8e3; fs = 6e3;
N = 1000;
n = 0:100000;
 xn = sin(2*pi*f0*n./fs+pi/2);
 playsound(xn,fs);
load('LPFlab1_1.mat');
% fc = pi/3/2/pi*fs;
% plotMagPhase(FTD((0:(length(lowPassFilter)-1)),lowPassFilter,N),N);
% freqz(lowPassFilter);
% 
% xnup(1:3:(3*length(xn)))= xn;
% xndown = xn(1:2:end);


%% Option 1
% A/D Converter
xn = sin(2*pi*f0*n./fs+pi/2);
%playsound(xn,fs);
% SRD
y1n = xn(1:3:end);
% % D/A Converter
 playsound(y1n,fs);
 playsound(y1n,fs/3);


%% Option 2
% A/D Converter
xn = sin(2*pi*f0*n./fs+pi/2);
% LPF
vn = conv(xn,lowPassFilter);
% SRD
y2n = vn(1:3:end);
% D/A Converter
%playsound(y2n,fs);
%playsound(y2n,fs/3);
figure(1);
plotMagPhase(FTD((0:(length(xn)-1)),xn,N),N);
figure(2);
plotMagPhase(FTD((0:(length(y1n)-1)),y1n,N),N);

%% Option 3
% A/D Converter
xn = sin(2*pi*f0*n./fs+pi/2);
% SRI
y3n(1:3:3*length(xn)) = xn;
% D/A Converter
%playsound(y3n,fs);
%playsound(y3n,fs*3);

%% Option 4
% A/D Converter
xn = sin(2*pi*f0*n./fs+pi/2);
% SRI
vn(1:3:3*length(xn)) = xn;
% LPF
y4n = conv(vn,lowPassFilter);
% D/A Converter
%playsound(y4n,fs);
%playsound(y4n,fs*3);


function [  ] = playsound(x,fs)
    sound(x,fs);
    pause(length(x)*1/fs+.5);
end