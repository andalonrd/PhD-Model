function [f] = LegQDs(v,x )

%Legendre second type function and its derivatives up the third. 


f(1,1) = LegQ(v,x);
f(2,1) = (v/(1-x^2))*(LegQ(v-1,x)-x*LegQ(v,x));
f(3,1) = (2*v*x/((1-x^2)^2))*LegQ(v-1,x)-(v*((1-v)*x^2+v+1)/((1-x^2)^2))*LegQ(v,x);
end

