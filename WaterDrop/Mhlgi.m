function [ h ] = Mhlgi(t,x)

h=(acos(x)/sin(acos(x)))^(1/2)*besseli(0,t*acos(x));