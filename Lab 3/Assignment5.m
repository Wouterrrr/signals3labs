clear all; close all;

%% Parameters
L = 10000;
ARorder = 6;

%% generate random data
% definition of filter
H.num = [1 0.56 0.81];
H.den = [1 -0.24 0.08 -0.37 0.52];

% definition of i
i = randn(10000,1);

% definition of x
x = filter(H.num, H.den, i);

%% Generate the true spectrum
%definition of power spectral
P.num = conv(H.num,H.num);
P.den = conv(H.den,H.den);

%setting information for plots
theta = 0:0.001:pi;
z = exp(1i.*theta);
n = 0:(length(P.den)-1);
Pden = P.den' .* z.^(-n');
Pz = 1./sum(Pden,1);

%actual plotting
plot(theta,20*log10(abs(Pz)))
hold on;

%% Blackman-Tukey' Method
% Generate the window function
% Window
alfa    = .54;
beta    = .46;
n = 0:(L-1);
Nw      = length(n);
rwtau       = [alfa - beta*cos(2*pi*n/(Nw - 1))'; zeros(length(x)-Nw,1)];

% calculate biased AC estimator
N = length(x);
rb.tau = (0:(N-1));
for tau = rb.tau
    rb.rb(tau+1,1) = 1/N .* sum(x(tau+1:end).*x(((tau+1):N) - tau));
end

% Calculate estimate
theta = 0:0.001:pi;
Pbt = sum(rwtau.*rb.rb.*exp(-1j.*rb.tau'.*theta),1);

% Plot the result
hold on;
% G = 16;
% Pbt = conv(1/G.*ones(1,G),abs(Pbt));
% Pbt = Pbt(1:end-(G-1));
plot(theta,20*log10(abs(Pbt)));


%% AR modeling method.
%find approximate AR values
A = aryule(x,ARorder);

% find the spectral power on basis of found AR values
P.num = 1;
P.den = conv(A,A);

% make the data ready for plottoing
theta = 0:0.001:pi;
z = exp(1i.*theta);
n = 0:(length(P.den)-1);
Pden = P.den' .* z.^(-n');
Pz = 1./sum(Pden,1);

% plot
hold on;
plot(theta,20*log10(abs(Pz)));
xlabel('\theta');
ylabel('P(e^{j\theta})');