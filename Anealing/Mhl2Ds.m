function [ f ] = Mhl2Ds(t,x)

%Legendre function and its derivatives up the third. 

fa = MhlDs(t,-1*x);
f(1,1) = fa(1,1);
f(2,1) = -fa(2,1);
f(3,1) = fa(3,1);
f(4,1) = -fa(4,1);

end

