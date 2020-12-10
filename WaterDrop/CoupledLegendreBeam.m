
function [ Aout, tel ] =CoupledLegendreBeam(X,alp,d,Tn,dmp)

%This function calculates the mode shapes for a given set of lenght 
%coordinates. FInally collects all the modal information into a single file

%Input Parameters:
% X = lenght Coordinates
% alp = alpha value (shear beam inf, flexural beam zero)
% d = delta (ratio of stiffness at top to the base)
% Tn = first Mode Period
% dmp = damping ratio

%Output Parameters
%Aout = File with the following strucutre:
% Fist Line, Periods of each mode
% Second Line, Damping ratio of each mode
% from there, the modal ordinates required at X coordinates

%tel = Running time of the procedure. 

tic;
warning('off','MATLAB:nearlySingularMatrix');

if d==0;

    d=0.0001;

end


T1 = importdata('T1.txt');

%Imports the Normalized periods computed with Mathematica)



dalp = T1(1,3)-T1(1,2);
dd = T1(3,1)-T1(2,1);
%differential ratio betwen sucesive alpha and damping 

T2 = importdata('T2.txt');
T3 = importdata('T3.txt');
T4 = importdata('T4.txt');
T5 = importdata('T5.txt');

%Number of alpha and damping values considered
Nalp = round(alp/dalp)+2;
Nd = round(d/dd)+2;


Ts(1) = T1(Nd,Nalp);
Ts(2) = T2(Nd,Nalp);
Ts(3) = T3(Nd,Nalp);
Ts(4) = T4(Nd,Nalp);
Ts(5) = T5(Nd,Nalp);

clear('T1','T2','T3','T4','T5');

%Imports the Top Normalized MPF values

MPF1 = importdata('MPF1.txt');
MPF2 = importdata('MPF2.txt');
MPF3 = importdata('MPF3.txt');
MPF4 = importdata('MPF4.txt');
MPF5 = importdata('MPF5.txt');

MPF(1) = MPF1(Nd,Nalp);
MPF(2) = MPF2(Nd,Nalp);
MPF(3) = MPF3(Nd,Nalp);
MPF(4) = MPF4(Nd,Nalp);
MPF(5) = MPF5(Nd,Nalp);

clear('MPF1','MPF2','MPF3','MPF4','MPF5');


Nx = length(X);
Aout = zeros(Nx+2,5);

if X(Nx) == 1
    
    X(Nx)=0.99999;

end

%if only one damping value is imposed, considers an uniform quantity for
%all modes
if length(dmp) == 1
    
    dmp =ones(5,1)*dmp;

end

for j=1:5

%computation of periods, obtained from the ratio of the normalized values imported   
Aout(1,j) = Tn*Ts(j)/Ts(1);    
Aout(2,j) = dmp(j);

[Pr A] = MCE(alp,d,Ts(j));

x =0.999999*sqrt(1-d);
    
Ps = real(LegP(Pr(1,1),x));
Qs = real(LegQ(Pr(1,1),x));
Mh1 = real(Mhlg(Pr(1,2),x,-0.5));
Mhl2 = real(Mhlg(Pr(1,2),-x,-0.5));
    
ModeL = A(1,1)*Ps(1,1)+A(2,1)*Qs(1,1)+A(3,1)*Mh1(1,1)+Mhl2(1,1);


    for i=1:Nx

        if X(i) == 1
           
            X(i) = 0.999999;
            
        end
        
        x = X(i)*sqrt(1-d);
    
        Ps = real(LegP(Pr(1,1),x));
        Qs = real(LegQ(Pr(1,1),x));
        Mh1 = real(Mhlg(Pr(1,2),x,-0.5));
        Mhl2 = real(Mhlg(Pr(1,2),-x,-0.5));
    
        Aout(i+2,j) = MPF(j)/ModeL*(A(1,1)*Ps(1,1)+A(2,1)*Qs(1,1)+A(3,1)*Mh1(1,1)+Mhl2(1,1));

    end
 
 tel = toc;   
    
end



