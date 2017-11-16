clear; close all; clc
%% Units
s   = 1;
Hz  = 1/s;
kHz = 10^3*Hz;
%% Load filter
load('LPF1.mat');
%% Parameters
f0 = 0.18*kHz;
fs = 6.0*kHz;
tstart  = 0;
tend    = 10;
%% 1) Create sampled sine 
x = zeros(3,tend*fs);
for i = 1:tend*fs;
    t       = (i-1)*1/fs + tstart;
    x(1,i)  = t;
    x(2,i)  = i - 1;
    x(3,i)  = sin(2*pi*f0*t + pi/2);
end
%% 1) Plotting
figure(1)
hold on
grid on
stem(x(2,1:50),x(3,1:50))
plot(fs*linspace(x(1,1),x(1,50),50*100), sin(2*pi*f0*linspace(x(1,1),x(1,50),50*100) + pi/2), 'k:')
legend('Sampled sinusoidal input','Underlying sine wave at t = n*fs')
xlabel('Index n')
ylabel('Sampled output x[n]')
%% 2) LPF
figure(2)
hold on
grid on
stem(LPF)
xlabel('Index n')
ylabel('Impulse response h[n]')

figure(3)
hold on
grid on
freqz(LPF)
subplot(2,1,1)
line([0.3*pi 0.3*pi],[-inf inf], 'Color', 'k')
%% 3) 
% y1
y1(1,:) = x(1,mod(x(2,:),3) == 0);                                          % Just copy time index of sample points
y1(2,:) = x(2,mod(x(2,:),3) == 0)/3;                                        % Indices are divided by 3
y1(3,:) = x(3,mod(x(2,:),3) == 0);                                          % Values of x[n] are copied once every 3 samples
% y2
y2(1,:) = x(1,mod(x(2,:),3) == 0);  
y2(2,:) = x(2,mod(x(2,:),3) == 0)/3;
output  = conv(LPF,x(3,:));                                                 % Output is x convolved with LPF
y2(3,:) = output(1,mod(x(2,:),3) == 0);                                         % y2 is output downsampled by factor 3    

