close all; clear all;

%% Assignment 1: Testing the cross-correlation function with energy signals
% Parameters
x1 = [1 2 3 2];
x2 = [3 1 2 1];
L = 5;

% Calculations
[l, Rxy] = crosscorxy(x1,x2,L);
[Rxy2,l2] = xcorr(x1,x2);

% Plots
stem(l,Rxy);
hold on;
stem(l2,Rxy2);
legend('self-created function', 'matlab function');
xlabel('l');
ylabel('$\hat{r}_{x,y}(l)$');

%c) The self created function is a normalized version of the
%cross-correlation. The matlab function is not normalized to the amount of
%samples, which increases the amplitude significantly. The matlab function
%does not correspond to the energy for l = 0


%% Assignment 2: Testing the cross-correlation function with random signals
clear all;
for k = 4
% creating samples
N = 10^k;
wn = normrnd(0,sqrt(2),1,N);

% Applying the filter
hn = [1 0.3 0.2];
xn = conv(wn(1:(end-2)),hn);

% Calculation Auto- and Crosscorrelation Estimates
L=10;
[lx,rx] = crosscorxy(xn,xn,L);
[lw,Rw] = crosscorxy(wn,wn,L);
[lxw,Rxw] = crosscorxy(xn,wn,L);
[lwx,Rwx] = crosscorxy(wn,xn,L);


% [Rx,lx] = xcorr(xn,xn);
% [Rw,lw] = xcorr(wn,wn);
% [Rxw,lxw] = xcorr(xn,wn);
% [Rwx,lwx] = xcorr(wn,xn);


% Plots
figure(k+1);
subplot(221);
stem(lx,rx);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{x}(l)');
subplot(222);
stem(lw,Rw);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{w}(l)');
subplot(223);
stem(lxw,Rxw);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{xw}(l)');
subplot(224);
stem(lwx,Rwx);
xlabel('l');
ylabel('r_{wx}(l)');
xlim([-L L]);
grid on;

end

% Question E
% it represents sigma_w^2, which is equal to 2 in this case, so it does
% correspond to the graph approximately, but the length is finite, so it is
% not a perfect match

% Bonus
% creating samples
un = [1 zeros(1,15)];

% Applying the filter
hn = [1 0.3 0.2];
xn = conv(un,hn);

%plot
figure;
stem(1:length(xn),xn);

%% Assignment 3: Experimenting with the cross-correlation function
clear all;
% Load data
load('exercise_correlation.mat');

% Calculate Correlation Functions
L=10;
[lx,rx] = crosscorxy(x,x,L);
[ly,Ry] = crosscorxy(y,y,L);
[lxy,Rxy] = crosscorxy(x,y,L);

% Plots
fig = figure;
subplot(131);
stem(lx,rx);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{x}(l)');
subplot(132);
stem(ly,Ry);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{y}(l)');
subplot(133);
stem(lxy,Rxy);
xlabel('l');
xlim([-L L]);
grid on;
ylabel('r_{xy}(l)');

% From the plots it can be concluded that rx has no repetitive components.
% ry does have repetitive components can be concluded from the second plot.
% The third plot shows that rxy has repetitive components and they seem to
% be causal only. 

%% Assignment 4: The channel
clear all; close all;
% a)
% definition of h
a = 2/5;
hn = [0 a 1 a];
n = 0:3;
% plot the frequency response
N = 1001;
plotMagPhase(FTD(n,hn,N),N);

% b)
figure;
zplane(hn);

% c) There is zero outside of the unit circle, so when you invert this, a
% pole outside of the unit circle will be there, which means that there
% will not be a stable system

%% Assignment 5: Ideal case and no channel noise, thus g[n] = 0
%define constants: note: A0 actually is equal to A0 = 5/2z^2
A0prime = 5/2; A1 = -1/3; A2 = 4/3; z1 = -1/2;
% define plot length
n = -8:8;
% calculate wn
wn = A0prime .* (A1.*(z1).^(n+1) + A2.*(1./z1).^(n+1));
% plotting
fig = figure;
plot(n,wn);

%testcase
xn = sin(pi/10.*n);
N = 10^5;
n = -N:N;
wn = A0prime .* (A1.*(z1).^(n+1) + A2.*(1./z1).^(n+1));
rn = conv(xn,hn);
yn = conv(rn,wn);
figure;stem(yn);
figure;stem(xn);

%% Assignment 6: FIR solution using optimal filter design
clear all; close all;
Dmax = 9;
% definition of h
a = 2/5;
hn = [0 a 1 a];
% definition of signals s[n] and g[n]
Varg = 0;
Vars = 1;
sigma_g = sqrt(Varg);
sigma_s = sqrt(Vars);
% definition of autocorrelations
L = 5;
l = -L:L;
Rg = zeros(1,length(l));
Rg(l==0) = sigma_g^2;
Rs = zeros(1,length(l));
Rs(l==0) = sigma_s^2;
% Calculations
rx = [0 0 Rg 0 0] + 1.32*[0 0 Rs 0 0] + 0.8*[0 Rs 0 0 0] + ...
     + 4/25*[Rs 0 0 0 0] + 0.8*[0 0 0 Rs 0] + 4/25*[0 0 0 0 Rs]; 
rx2 = conv(hn, fliplr(hn));
l = -(L+2):(L+2);
% plotting
stem(l,rx);

%setting up the matrix
S = 11;
Rx = zeros(S,S);
n = diag(ones(1,S));
for N = 1:S
    for M = 1:S
        if M==N
            Rx(N,M) = rx(l==0);
        elseif M == N+1 || N == M+1
            Rx(N,M) = rx(l==1);
        elseif M == N+2 || N == M+2
            Rx(N,M) = rx(l==2);
        end
            
    end
end

% calculation of rdx
rs0 = sigma_s^2;
a = 2/5;
ldx = -2:10;
n = -100:100;
delta = zeros(1,length(n));
delta(n==0) = 1;
x = 0:(S-1);
rdx = zeros(Dmax+1,S);
for D = 0:Dmax
    for l = x
        rdx(D+1,x==l) = a*delta(n==l-D+1) + delta(n==l-D+2) + a*delta(n==l-D+3);
        
    end  
end

%calculation of optimum filter
for D = 0:Dmax
   hopt(D+1,:) = (inv(Rx)*((rdx(D+1,:))'))'; 
end

% implementation of system
a = 2/5;
Varg = 0.5;
N = 10^4;
hn = [0 a 1 a];
for D = 0:Dmax
   wDn = conv(hn,hopt(D+1,:));
   subplot(4,5,D+1);
   stem(wDn);
   hold on;
end

for D = 0:Dmax
    MSE(D+1) = Vars - rdx(D+1,:)*hopt(D+1,:)';
end
fig = figure;
stem(0:length(MSE)-1,20*log(MSE));
xlabel('D');
ylabel('MSE [dB]');
grid on;
saveas(fig,'Assignment6d.png');


fig = figure;
n = -8:8;
w_id = 5/2*2/3*( (-1/2).^(n+2).*Step(n+2) + (-2).^(n+2).*Step(-n -3) );
for D = 0:Dmax
   subplot(4,3,D+1);
   stem(0:S-1,hopt(D+1,:));
   %hold on;
   %stem(n+D, w_id);
   title(['D = ' num2str(D)]);
   ylabel('');
   xlabel('n');
   grid on;
end
[ax2,h2]=suplabel('w_D[n]','y');
saveas(fig,'Lab2Assigment6b1.png');

fig=figure;
for D = 0:Dmax
   subplot(4,3,D+1);
   stem(0:S+2,conv(hopt(D+1,:),hn));
   title(['D = ' num2str(D)]);
   xlabel('n');
   grid on;
end
[ax2,h2]=suplabel('g_D[n] = h[n] * w_D[n]','y');

saveas(fig,'Lab2Assigment6b2.png');
function [ u ] = Step(n)
    u = double(n >= 0);
end