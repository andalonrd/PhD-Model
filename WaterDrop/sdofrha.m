% Function to compute the response history analysis using recursive method
% Pages 167-171 of Chopra's book - 2nd edition

% Original reference:
% Nigam and Jennings (BSSA, v. 59, 909-922, 1969)

% Required data
% accel = vector containing the acceleration time history
% per = period of vibration of the SDOF system
% xi = damping ratio of the SDOF system
% d0 = initial displacement
% v0 = initial velocity

% Data returned by the function
% rd = vector containing the relative displacement time history
% rv = vector containing the relative velocity time history
% ra = vector containing the relative acceleration time history
% aa = vector containing the absolute acceleration time history

function [rd,rv,ra,aa,f]=sdofrha(accel,per,xi,dt,d0,v0)

% initializes vectors where response parameters will be stored
rd=zeros(size(accel));
rv=zeros(size(accel));
ra=zeros(size(accel));
aa=zeros(size(accel));

% computes some parameters of the SDOF system
wn = 2*pi/per;
wd = wn*(1-xi*xi)^0.5;
k = wn*wn*1.0;

% computes coefficient for recirrence equations
A = (exp(-xi*wn*dt))*((xi/sqrt(1-xi*xi))*sin(wd*dt)+cos(wd*dt));
B = (exp(-xi*wn*dt))*((1/wd)*sin(wd*dt));
C = (1/k)*((2*xi/wn/dt)+(exp(-xi*wn*dt))*(((1-2*xi*xi)/wd/dt-xi/sqrt(1-xi*xi))*sin(wd*dt)-(1+2*xi/wn/dt)*cos(wd*dt)));
D = (1/k)*(1-(2*xi/wn/dt)+(exp(-xi*wn*dt))*(((2*xi*xi-1)/wd/dt)*sin(wd*dt)+(2*xi/wn/dt)*cos(wd*dt)));
Ap = -(exp(-xi*wn*dt))*((wn/sqrt(1-xi*xi))*sin(wd*dt));
Bp = (exp(-xi*wn*dt))*(cos(wd*dt)-(xi/sqrt(1-xi*xi))*sin(wd*dt));
Cp = (1/k)*((-1/dt)+(exp(-xi*wn*dt))*(((wn/sqrt(1-xi*xi))+xi/dt/sqrt(1-xi*xi))*sin(wd*dt)+(1/dt)*(cos(wd*dt))));
Dp = (1/k/dt)*(1-((exp(-xi*wn*dt))*(xi/sqrt(1-xi*xi)*sin(wd*dt)+cos(wd*dt))));

% response history analysis
rd(1)=d0;
rv(1)=v0;
p(1)=-accel(1);
ra(1)=accel(1);

for i=2:(length(accel)-1);
    p(i)=-accel(i);
    rd(i)=A*rd(i-1)+B*rv(i-1)+C*p(i-1)+D*p(i);
    rv(i)=Ap*rd(i-1)+Bp*rv(i-1)+Cp*p(i-1)+Dp*p(i);
    ra(i)=p(i)-2*xi*wn*rv(i)-k*rd(i);
    aa(i)=accel(i)+ra(i);
end

f=k*rd;

return