function [ output_args ] = plotMagPhase( freqResponse,N )
%PLOTMAGPHASE Summary of this function goes here
%   Detailed explanation goes here
    theta = -pi:2*pi/(N-1):pi;
    m = abs(freqResponse);                               % Magnitude
    p = unwrap(angle(freqResponse));                     % Phase
    subplot(2,1,1);
    %hold on;
    plot(theta,m)
    title('Magnitude');
    ylabel('|H_{causal}(e^{j\theta}| [-]');
    xlabel('\theta [rad/s]');
    ax = gca;
    ax.XTick = [-pi,-pi/2,0,pi/2,pi];
    grid on;

    subplot(2,1,2)
    plot(theta,p)
    title('Phase')
    ylabel('\phi_{H_{causal}}(e^{j\theta}) [rad]');
    xlabel('\theta [rad/s]');
    ax = gca;
    ax.XTick = [-pi,-pi/2,0,pi/2,pi];
    grid on;

end

