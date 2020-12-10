function [ h ] =d2Mhlgi(t,x)

h = real((d1Mhlgi(t,1.001*x)-d1Mhlgi(t,0.999*x))/((1.001-0.999)*x));

end

