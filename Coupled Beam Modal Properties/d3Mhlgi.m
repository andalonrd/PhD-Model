function [ h ] = d3Mhlgi( t,x )

h = real((d2Mhlgi(t,1.001*x)-d2Mhlgi(t,0.999*x))/((1.001-0.999)*x));


end

