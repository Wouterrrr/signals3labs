function [l, Rxy] = crosscorxy(x1,x2,L)
    
    %make sure data has right dimensions
    if size(x1,2) == 1
        x1 = x1';
    end
    if size(x2,2) == 1
        x2 = x2';
    end

    l =-L:L;
    zeroPadLength = max(length(x1),L);
    X1 = [zeros(1,zeroPadLength) x1 zeros(1,zeroPadLength)];
    X2 = [zeros(1,zeroPadLength) x2 zeros(1,zeroPadLength)];

    if length(x1) ~= length(x2)
       error('Length of input signals should be equal.'); 
    end
    
    for i = 1:length(l)
        li = l(i);
        Rxy(i) = 0;
        for n = 1:length(x1)
           Rxy(i) = Rxy(i) + X1(n+zeroPadLength)*conj(X2(n+zeroPadLength-li));
        end
    end
    Rxy = Rxy/length(x1);
end

