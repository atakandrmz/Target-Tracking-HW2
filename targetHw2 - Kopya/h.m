function y = h(x,y, measNoise)
   % g = [sin(x)/(10.*sqrt(x.^2+y.^2)); cos(y)/(10.*sqrt(x.^2+y.^2))]
    y = [x y] + measNoise
end