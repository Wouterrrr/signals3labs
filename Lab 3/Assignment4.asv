clear all; close all;
%% parameters
ARorder = 100;

%% generate random data
% definition of filter
H.num = 1;
H.den = [1 -1.5 0.64];

% definition of i
i = randn(1000,1);

% definition of x
x = filter(H.num, H.den, i);

%% Actual Spectrum
%definition of power spectral
P.num = 1;
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

%% find AR parameters
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
legend('Correct spectrum', 'Spectrum estimated with n = 2');

%% Blackman-Tukey' Method
L = 10;
% Generate the window function
rwtau = [ones(L,1); zeros(length(x)-L,1)];

% calculate biased AC estimator
N = length(x);
rb.tau = (0:(N-1));
for tau = rb.tau
    rb.rb(tau+1,1) = 1/N.* sum(x(tau+1:end).*x(((tau+1):N) - tau));
end

% Calculate estimate
theta = 0:0.001:pi;
Pbt = sum(rwtau.*rb.rb.*exp(-1j.*rb.tau'.*theta),1);

% Plot the result
hold on;
plot(theta,20*log10(abs(Pbt)));