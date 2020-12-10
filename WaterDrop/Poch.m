function [ h ] = Poch( A,B )
%Pochhammer Symbol (in numeric terms and complex values)
% 
%Input parameters:
% A First Parameter, real or complex
% B Term of the Pochhammer simbol (integer number)
% h Pochhammer 

h = gammac(A+B)/(gammac(A));

end

