function [Vt X] = MCE(alp,d,Te )

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

wt = 2*pi()/Te;
alpe = alp/(1-d)^0.5;
wte = wt/(1-d);
bo = (2-alpe^2+((alpe^2-2)^2+4*wte^2)^0.5)/2;
ao = (2-alpe^2-((alpe^2-2)^2+4*wte^2)^0.5)/2;
ld = (abs(ao)-0.25)^0.5;
v  = -0.5+(1+4*bo)^0.5/2;

A = LegPDs(v,0);

if d<0.0001 
B = LegPDs(v,0.9999);
else
B = LegPDs(v,(1-d)^0.5);
end

C = LegQDs(v,0.0001);

if d<0.0001 
D = LegQDs(v,0.9999);
else
D = LegQDs(v,(1-d)^0.5);
end

E = MhlDs(ld,0.0001);

if d <0.0001 
F = MhlDs(ld,0.9999);
else
F = MhlDs(ld,(1-d)^0.5);
end

G(1,1) = E(1,1);
G(2,1) = -E(2,1);

if d <0.0001 
H = Mhl2Ds(ld,0.9999);
else 
H = Mhl2Ds(ld,(1-d)^0.5);    
end

M(1,1) = A(1,1);
M(1,2) = C(1,1);
M(1,3) = E(1,1);
M(1,4) = G(1,1);

M(2,1) = A(2,1);
M(2,2) = C(2,1);
M(2,3) = E(2,1);
M(2,4) = G(2,1);

M(3,1) = B(3,1);
M(3,2) = D(3,1);
M(3,3) = F(3,1);
M(3,4) = H(3,1);

A = M(1:3,1:3);
B = -M(1:3,4);
X = A\B;
Vt = [v ld];
end

