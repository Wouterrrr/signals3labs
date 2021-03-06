clear all; close all;

%% Original
%parameters
fs = 10e6;
sig = 3/fs;
N = 10;
M = 8;
k = 3;
alpha = 0.0001;

%definition of signal that needs to be shifted
n = -N:N;
xn = gaussmf(n*1/fs,[sig 0]);

%Definition Low Pass Filter
nh = -100000:100000;
LPF = M./M.*sinc(nh./M);

tic;
%upsampling
vn(1:M:(length(xn)*M)) = xn;
n = -N:1/M:N;

%delay
vn = [zeros(1,k) vn];
hold on;

%Applying Low Pass Filter on original signal
vn = conv(vn,LPF);

%shift signal back to the left as compensation for the causal Low Pass
%Filter, which introduced a delay of a few samples.
vn = vn((length(LPF)+1)/2:end);

%Downsampling
yn = vn(1:M:end);
originalMethodTime = toc;

%plots
fig = figure('NumberTitle', 'off', 'Name', 'Assignment 7: Original Method');
nplot = -N:(length(yn)-N-1);
stem(0:(length(yn)-1),yn);
hold on;
stem(0:(length(xn)-1),xn);
n = linspace(0,length(xn)-1,1000);
plot(n,gaussmf((n-10)*1/fs,[sig 0]));
plot(n,gaussmf((n-10-3/8)*1/fs,[sig 0]))
xlim([0 30]);
xlabel('n');
legend('shifted sampled signal', 'sampled signal', 'non-sampled signal', 'shifted non-sampled signal');
grid on;
saveas(fig,'Assignment7originalmethod.png');

%% Second option: PolyPhase
%definition signal that has to be shifted
nx = -10:10;
x = gaussmf(nx*1/fs,[sig 0]);

%definition of Low Pass Filter
nh = -100000:100000;
LPF = M./M.*sinc(nh./M);

%calculation Mth polyPhase components
polyPhase = LPF(1+M-k:M:end);

tic;
%Calculation of shifted version of signal
y = conv(polyPhase,x);

%because the filter is linear phase and causal, it causes a shift of a
%certain amount of samples to the right. The compensation for this happens
%here.
y = y((length(polyPhase)+1)/2-1:end);
polyPhaseMethodTime = toc;

%plots
fig = figure('NumberTitle', 'off', 'Name', 'Assignment 7: PolyPhase method');
stem(0:(length(y)-1),y);
hold on;
stem(0:(length(xn)-1),xn);
n = linspace(0,length(xn)-1,1000);
plot(n,gaussmf((n-10)*1/fs,[sig 0]));
plot(n,gaussmf((n-10-3/8)*1/fs,[sig 0]))
xlim([0 30]);
xlabel('n');
legend('shifted sampled signal', 'sampled signal', 'non-sampled signal', 'shifted non-sampled signal');
grid on;
saveas(fig,'Assignment7polyphase.png');
