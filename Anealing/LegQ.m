
function [ h ] = LegQ(v,x)
%
% Legendre function of the second kind and degree v

% Parameters
% v = Degre of legendre function
% x = ordinate

% Output
% h = Function Output

h = pi()/2*(LegP(v,x)*cos(v*pi())-LegP(v,-1*x))/sin(v*pi());
end

