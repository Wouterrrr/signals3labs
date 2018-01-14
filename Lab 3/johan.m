clear; close all; clc
%% Assignment 1
%% 1
% System function
a = 1.5;
b = 0.64;
Hz.num = [1];
Hz.den = [1, -a, b];
% PSD
fs = 1;
f = 0:0.001:.5;
% such that f0 = f/fs = 0:0.001:.5
Px.num = [0, 0, 1];
Px.den = [b, -a*(b+1), (1+a^2+b^2), -a*(b+1), b];
Po = freqz(Px.num, Px.den, f, fs);

figure(1)
hold on
grid on
hpo = plot(f,20*log10(abs(Po)));
xlabel('Normalized frequency [-]')
ylabel('Magnitude/Power [dB]')

legendstr1{1,1} = 'True spectrum P_x(e^{j\theta})';
plotmatrix1(1,1) = hpo;

figure(2)
hold on
grid on
hpo = plot(f,20*log10(abs(Po)));
xlabel('Normalized frequency [-]')
ylabel('Magnitude/Power [dB]')

legendstr2{1,1} = 'True spectrum P_x(e^{j\theta})';
plotmatrix2(1,1) = hpo;
%% 2
% N = 1024;
N = 2^10;
L = 50;
x = randn(N+L,1);
x = x(L+1:end);
x = filter(Hz.num, Hz.den, x);
padfactor = 0;
X    = fft([x; zeros(padfactor*size(x,1), size(x,2))]);
P    = 1/N*X.*conj(X);
range= 0:1/(N*(padfactor+1)):(N-1)/(2*N);
figure(1)
hp = plot(range, 20*log10(P(1:length(range))));
%hp = plot(-(N-1)/(2*N):1/N:(N-1)/(2*N),20*log10(P));

legendstr1{1,length(legendstr1)+1} = 'Periodogram P';
plotmatrix1(1,length(plotmatrix1)+1) = hp;

%% 3
ii = 1;
for K = [4, 16]
    M(ii) = N/K;
    range = 0:K/((padfactor+1)*N):(N-1)/(2*N);
    x_segmented{ii} = reshape(x,M(ii),K);
    X_segmented{ii} = fft([x_segmented{ii}; zeros(padfactor*size(x_segmented{ii},1), size(x_segmented{ii},2))]);
    P_segmented{ii} = 1/M(ii)*X_segmented{ii}.*conj(X_segmented{ii});
    PB{ii}          = mean(P_segmented{ii},2);
    figure(2)
    hb = plot(range, 20*log10(PB{ii}(1:length(range))));
    legendstr2{length(legendstr2)+1} = ['Bartlett K = ', num2str(K)];
    plotmatrix2(length(plotmatrix2)+1) = hb;
    ii = ii + 1;  
end
%% Assignment 2
%% 4
load('Window2_4.mat'); %contains Hamming window of length 64
PW=pwelch(x,window_1,[ ],2*pi*f)*2*pi;
figure(2)
hpw = plot(f,20*log10(PW));
legendstr2{length(legendstr2)+1} = 'Welch';
plotmatrix2(length(plotmatrix2)+1) = hpw;
%% Assignment 3
%% 1: See function
%% 2
M = length(x)/4;
NFFT = length(x);
[PBT,fgrid] = btmethod(x,M,NFFT);
figure(2)
hbt = plot(fgrid,20*log10(abs(PBT)));
legendstr2{length(legendstr2)+1} = 'Blackman-Turkey';
plotmatrix2(length(plotmatrix2)+1) = hbt;

figure(3)
grid on
hold on
hpo = plot(f,20*log10(abs(Po)));
legendstr3{1} = 'True spectrum P_x(e^{j\theta})';
plotmatrix3(1) = hpo;
xlabel('Normalized frequency [-]')
ylabel('Magnitude/Power [dB]')
hbt = plot(fgrid,20*log10(abs(PBT)));
legendstr3{length(legendstr3)+1} = 'Blackman-Turkey';
plotmatrix3(length(plotmatrix3)+1) = hbt;

%% Update legend
figure(1)
legend(plotmatrix1, legendstr1, 'Location', 'NorthEast')
figure(2)
legend(plotmatrix2, legendstr2, 'Location', 'NorthEast')
figure(3)
legend(plotmatrix3, legendstr3, 'Location', 'NorthEast')
%% Functions
function [PBT,fgrid] = btmethod(x,M,NFFT)
%% btmethod calcalutes the spectral estimate PBT and the corresponding 
% frequency grid fgrid. Takes data x, upper number of lags M and the FFT
% length NFFT. Window is Hamming with length 2M + 1, so M samples to the
% left and M to the right.
% For step a) use hamming window of length 2M+1, and in step b) you can use the command
% rx=xcorr(x,M,'biased').
% Check input
assert(M>=1, 'M should be bigger than zero as it is a measure for the one sided window length!');
assert(NFFT >= 2*M+1, ['NFFT should at least be equal to the window length (2*M+1) to ' ...
        'at least transform just the autocorrelation. NFFT > M means' ...
        'padding with NFFT - M zeros.'])
% Window
alfa    = .54;
beta    = .46;
n       = 0:2*M;
Nw      = length(n);
assert(Nw == 2*M+1)
w       = alfa - beta*cos(2*pi*n/(Nw - 1))';
% Biased autocorrelation estimate
rx      = xcorr(x,M,'biased');
% Spectral estimate
PBT     = fft([rx.*w; zeros(NFFT-Nw,1)]);
% Corresponding frequency range
Nrx     = length(rx);
padfactor   = length(PBT)/Nrx - 1;
fgrid   = 0:1/(Nrx*(padfactor+1)):(Nrx-1)/(2*Nrx);
% Only for f0 = 0 to 0.5, so theta = 0 to pi and leave out the symmetric
% part f0 = -0.5 to 0 (theta = -pi to 0)
PBT     = PBT(1:length(fgrid));
end