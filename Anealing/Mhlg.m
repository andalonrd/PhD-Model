
function [ h ] = Mhlg(t,x,v)

% Conical Mehler function of order t

% Parameters
% x = ordinate where the function is being evaluated
% l= l parameter of Mehler Function
% h = Mehler funtion expressed in terms of Gauss Hypergeometric functions
% NOTE: v parameter must be -1/2! for implementation reasons, it is in the function set as variable, 
% but the function must be called always with v=-1/2! 

if t<15

h = LegP(v+sqrt(-1)*t,x);

else
    
h = real(Mhlgi(t,x));    

end

