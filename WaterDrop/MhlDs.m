function [ f ] = MhlDs(t,x)

%Legendre function and its derivatives up the third. 

if t<15

v=-0.5+t*1i; 

f = Mhlg(t,x,-0.5);
d1 = (v/(1-x^2))*(Mhlg(t,x,-0.5-1)-x*Mhlg(t,x,-0.5));
d2 = (2*v*x/((1-x^2)^2))*Mhlg(t,x,-0.5-1)-(v*((1-v)*x^2+v+1)/((1-x^2)^2))*Mhlg(t,x,-0.5);
d3 = (2*v*(v-1)*x/((1-x^2)^3))*Mhlg(t,x,-0.5-2) +((8-2*v)*v*x^2/((1-x^2)^3)+(2-2*v)*v/((1-x^2)^2)-v*(v-1)*((2-v)*x^2+v)/((1-x^2)^3))*Mhlg(t,x,-0.5-1)-x*((2-2*v)*v/((1-x^2)^2)+(4-v)*v*((1-v)*x^2+v+1)/((1-x^2)^3))*Mhlg(t,x,-0.5);

f(1,1) = real(f);
f(2,1) = real(d1);
f(3,1) = real(d2);
f(4,1) = real(d3);

else

f(1,1) = Mhlgi(t,x);
f(2,1) = d1Mhlgi(t,x);
f(3,1) = d2Mhlgi(t,x);
f(4,1) = d3Mhlgi(t,x);
end

