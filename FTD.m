function [ X ] = FTD( nx, x, N )
%FTD Summary of this function goes here
%   Detailed explanation goes here
    theta = -pi:2*pi/(N-1):pi;
    X = zeros(1,N);
    for i = 1:length(x)
       X = X + x(i)*exp(-1j.*theta.*nx(i));
    end

end

