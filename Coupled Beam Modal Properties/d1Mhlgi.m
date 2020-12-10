function [ h ] = d1Mhlgi( t,x )

h = real((Mhlgi(t,1.001*x)-Mhlgi(t,0.999*x))/((1.001-0.999)*x));

end

