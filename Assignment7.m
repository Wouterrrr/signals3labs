clear all; close all;

%% Original
%parameters
fs = 10e6;
sig = 3/fs;
N = 10;
M = 8;
k = 3;
alpha = 0.0001;

%pulse
n = -N:N;
xn = gaussmf(n*1/fs,[sig 0]);
figure;
stem(0:(length(xn)-1),xn);

%upsampling
vn(1:M:(length(xn)*M)) = xn;
n = -N:1/M:N;
figure;
stem(0:(length(vn)-1),vn);

%delay
vn = [zeros(1,k) vn];
hold on;
n = [n (N+1):(N+3)];
stem(0:(length(vn)-1),vn);

%Low Pass Filter
LPF = [0 0];
i = 0;
while mod(length(LPF),2) == 0
    LPF = fdesign.lowpass('Fp,Fst,Ap,Ast',0.99/M,1.01/M,0.1,0.001+alpha*i,'linear');
    Hd = design(LPF,'equiripple');
    LPF = Hd.Numerator;
    i = i+1;
end
vn = 8*conv(vn,LPF);
vn = vn((length(LPF)+1)/2:end);

figure;
stem(0:(length(vn)-1),vn);

%Downsampling
yn = vn(1:M:end);

%plotting
figure;
nplot = -N:(length(yn)-N-1);
stem(0:(length(yn)-1),yn);
hold on;
n = -N:N;
stem(0:(length(xn)-1),xn);
n = linspace(0,length(xn)-1,1000);
plot(n,gaussmf((n-10)*1/fs,[sig 0]));
plot(n,gaussmf((n-10-3/8)*1/fs,[sig 0]))

close all;
%% Second option
% polyPhase = LPF(1+M-k:M:end);
% vn = [0 xn];
% vn = 8*conv(vn,polyPhase);
% yn = vn((length(polyPhase)+1)/2+1:end);
% figure;
% stem(yn);
% 
% %plotting
% figure(6);
% hold off;
% stem(0:(length(yn)-1),yn);
% hold on;
% n = -N:N;
% stem(0:(length(xn)-1),xn);
% n = linspace(0,length(xn)-1,1000);
% plot(n,gaussmf((n-10)*1/fs,[sig 0]));
% plot(n,gaussmf((n-10-3/8)*1/fs,[sig 0]))

% clearvars -except M k;
nh = -100000:100000;
LPF = M./M.*sinc(nh./M);
%polyPhase = LPF(mod(-nh+5,8)==0);
polyPhase = LPF(1+M-k:M:end);
%npoly = (nh(mod(-nh+5,8)==0)-5)/8;
nx = -10:10;
x = gaussmf(nx*1/fs,[sig 0]);
y = conv(polyPhase,x);
y = y((length(polyPhase)+1)/2-1:end);
%[ny,y] = convcool(nx,x,npoly,polyPhase);
close all;
%plotting
stem(0:(length(y)-1),y);
hold on;
stem(0:(length(xn)-1),xn);
n = linspace(0,length(xn)-1,1000);
plot(n,gaussmf((n-10)*1/fs,[sig 0]));
plot(n,gaussmf((n-10-3/8)*1/fs,[sig 0]))
xlim([0 30]);

function [ ny, y ] = convcool(nx, x, nh, h  )
    y = conv(x,h);
    ny = linspace(nx(1) + nh(1), length(y)+nx(1) + nh(1) - 1,length(y));
end