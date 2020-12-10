function [ f ] = LegPDs(v,x )

%Legendre function and its derivatives up the third. 


f(1,1) = LegP(v,x);
f(2,1) = (v/(1-x^2))*(LegP(v-1,x)-x*LegP(v,x));
f(3,1) = (2*v*x/((1-x^2)^2))*LegP(v-1,x)-(v*((1-v)*x^2+v+1)/((1-x^2)^2))*LegP(v,x);

end

