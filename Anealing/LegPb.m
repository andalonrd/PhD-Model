

function [ h ] = LegPb(v,x )

% Legendre function of degree v (zero order)

% Parameters
% x = ordinate where the function is being evaluated
% v = degree of legendre function

% Output
% h = Legendre funtion evaluated for 75 terms

h = 0;

% Variables corresponding to each one of the terms of equation (4.6.1) in Zhang and Jin (1996)

if x<0 

    for k=0:75
       
        h =h+ Poch(-v,k)*Poch(v+1,k)/(factorial(k)^2)*(0.5+x/2)^k*(psic(k-v)+psic(k+v+1)-2*psic(k+1)+log(0.5+x/2));
        
    end
 h=sin(v*pi())/pi()*h;   

else

    for k=0:75
    
    h = h + Poch(-v,k)*Poch(v+1,k)/(factorial(k)^2)*(0.5-x/2)^k;
    
    end


end
end

