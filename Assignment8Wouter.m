clear all; close all;

%% Parameters
M = 10;
fs = 44100;
A=1;
N = 1e5;
n= 0:1e4;
[xn,fs] = audioread('beethoven5.wav');
xn = xn';
xn = cos(pi./2.*n);
WM = exp(-1j*2*pi/M);

%% Bandpass filter definitions
n = -N:N;
LPF = A.*sinc(n./M);

%bandpass filters
for g = 1:M
   BPF(g,:) = LPF.*cos(2*pi/M.*(g-1).*n); 
    
end

%% polyPhase
k = (1:M)';
figure;
hold on;
for i = k'
    polyPhase(i,:) = LPF(i:M:floor(length(LPF)/M)*M);
    stem(n(i:M:floor(length(LPF))),LPF(i:M:floor(length(LPF))));
end

%% apply polyphases to signal
lmax = length(polyPhase(1,:))+M-1+length(xn);
for g = k'
   vn = [zeros(1,(g-1)) xn];
   vn = vn(1:M:end);
   yn = conv(vn,polyPhase(g,:));
   un(g,:) = [yn zeros(1,lmax-length(yn))];
end

%% Matrix definition
i = 0:(M-1);H
p = i';
F_M = WM.^(i.*p);

%Matrix processing
vn = conj(F_M)*un;

%plotting
figure(3);
N = 1000;
plotMagPhase(FTD(0:length(un(1,:)),real(un(1,:)),N),N);figure(4);
plotMagPhase(FTD(0:length(un(2,:)),real(un(2,:)),N),N);figure(5);
plotMagPhase(FTD(0:length(un(3,:)),real(un(3,:)),N),N);

%% Processing / Equalizer
coeff = ones(1,M);
P = diag(coeff);
vn = P*vn;

%Matrix processing
vn = F_M * vn;

%Filter G
for g = k'
   yn = conv(vn(g,:),polyPhase(g,:));
end

%upsampling
for g = k'
   yn2(g,1:M:length(yn)*M) = yn(g);
end


%delays
yn = 0;
for g = k'
   yn = yn + [zeros(1,M-g) yn2(g,:) zeros(1,g)];
end

%plots
figure;
N = 1000;
plotMagPhase(FTD(0:length(xn),real(xn),N),N);
figure;
plotMagPhase(FTD(0:length(yn),real(yn),N),N);



