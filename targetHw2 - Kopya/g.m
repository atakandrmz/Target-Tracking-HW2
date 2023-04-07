function g = g(x,xv, y, yv, acc)
   % g = [sin(x)/(10.*sqrt(x.^2+y.^2)); cos(y)/(10.*sqrt(x.^2+y.^2))]
    g = [1 0.1 0 0;0 1 0 0; 0 0 1 0.1;0 0 0 1]*[x;xv;y;yv]+[0.01/2;0.1;0.01/2;0.1]*acc
end