function gBar = g(x,xv, y, yv, acc)
   % g = [sin(x)/(10.*sqrt(x.^2+y.^2)); cos(y)/(10.*sqrt(x.^2+y.^2))]
    gBar = [1 0.1 0 0;0 1 0 0; 0 0 1 0.1;0 0 0 1]
end