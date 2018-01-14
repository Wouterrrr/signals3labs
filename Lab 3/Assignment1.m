clear; close all;
figure(1); hold on;

%% a
f = 0:0.001:0.5; fs = 1; 
Px.num = [0 0 1];
Px.den = [0.64 -2.46 3.6596 -2.46 0.64];
z = exp(1i.*2.*pi.*f);
H = 1./(1 - 1.5.*z.^-1 + 0.64*z.^-2);
P = H.*conj(H);
plot(f, 20.*log10(abs(P)));
% 
% h = freqz(Px.num, Px.den, f,fs);
% xlabel('Normalized Frequency [-]');
% subplot(212);
% xlabel('Normalized Frequency [-]');

%% b
clear all;
N = 1024; L = 50;
H.num = 1;
H.den = [1 -1.5 0.64];
x = randn(N+L,1);
y = filter(H.num, H.den, x);
y = y((L+1):end,1);
n = 0:(length(y)-1);
Nplot = 10000;
P = FTD(n,y,Nplot);
P = P.*conj(P)/N;
plotMagPhase(P,Nplot);

%% c
% create segments
K = 1;
y = reshape(y',[K,N/K]);
n = 0:(length(y)-1);

for l = 1:K
    P = FTD(n,y(l,:),Nplot);
    PB(l,:) = P.*conj(P)/N;    
end
PB = mean(PB, 1);
plotMagPhase(PB, Nplot);






function [ X ] = FTD( nx, x, N )
%FTD Summary of this function goes here
%   Detailed explanation goes here
    f = 0:1/N:(N-1)/(2*N);
    X = zeros(1,N/2);
    for i = 1:length(x)
       X = X + x(i)*exp(-1j.*2*pi*f.*nx(i));
    end

end

function [ output_args ] = plotMagPhase( freqResponse,N )
%PLOTMAGPHASE Summary of this function goes here
%   Detailed explanation goes here
    f = 0:1/N:(N-1)/(2*N);
    m = 20*log10(abs(freqResponse));                               % Magnitude
    p = unwrap(angle(freqResponse));
    %hold on;
    plot(f,m)
    title('Magnitude');
    ylabel('|H_{causal}(e^{j\theta}| [dB]');
    xlabel('\theta [rad/s]');
    ax = gca;
    ax.XTick = [-pi,-pi/2,0,pi/2,pi];
    grid on;

%     subplot(2,1,2)
%     plot(f,p)
%     title('Phase')
%     ylabel('\phi_{H_{causal}}(e^{j\theta}) [rad]');
%     xlabel('\theta [rad/s]');
%     ax = gca;
%     ax.XTick = [-pi,-pi/2,0,pi/2,pi];
%     grid on;

end

