function [ Log ] = Anealing(T0,dmp0,alp0,d0)

[BaseMotion, Path] = uigetfile;
warning('off','MATLAB:nearlySingularMatrix');
filename = strcat(Path,BaseMotion);

A=importdata('T1.txt');
[Nds,NAlps]=size(A); 
dmax = A(Nds,1);
alpmax = A(1,NAlps);
dmin = A(2,1);
alpmin = A(1,2);
dd = A(3,1)-A(2,1);
dalp = A(1,3)-A(1,2);
clear('A');

InFile = importdata(filename);
[Nth Nx] = size(InFile);
X = InFile(2,3:Nx);
Ga = InFile(3:Nth,2);
AcTh = InFile(3:Nth,3:Nx);
Nx = Nx-2;

filename = strcat(Path,'d',BaseMotion);
InDisplacements = importdata(filename);
dThr = InDisplacements(3:Nth,:);
clear ('InDisplacements');
dTh = zeros(Nth-2,Nx);

for i=3:Nx+2;

    dTh(:,i-2)=dThr(:,i)-dThr(:,2);

end

dt = InFile(4,1)-InFile(3,1);
clear ('InFile');

Tmax = 6;
dT = 0.025;

dmpmax = 0.1;
ddmp = 0.0025;

dds(1,1) = dT;
dds(1,2) = ddmp;
dds(1,3) = dalp;
dds(1,4) = dd;

Mxs(1,1) = Tmax;
Mxs(1,2) = dmpmax;
Mxs(1,3) = alpmax;
Mxs(1,4) = dmax;

Mns(1,1) = 0;
Mns(1,2) = 0;
Mns(1,3) = alpmin;
Mns(1,4) = dmin;

B0(1,1) = T0;
B0(1,2) = dmp0;
B0(1,3) = alp0;
B0(1,4) = d0;
[Tms,dmp, Gphis,~] = CoupledParabolicBeam(X,B0(1,3),B0(1,4),B0(1,1),B0(1,2));
B0(1,5) = OneID(Ga,AcTh,dTh,dt,Tms,dmp,Gphis,0);
B0(1,6) = 0;

B1 = zeros(1,6);
Log = B0;

jmax = 2500;
ExitC = B0(1,5)/100;
LogN = 1;
Cmax = 3;
NBmax = 4*Cmax;
PostAnnealing = 0.05;
for j=1:jmax
    
    if rand(1) < PostAnnealing
        
        C = Cmax*(jmax-j)/jmax;
    
    else
            
        C= Cmax*rand(1);
        
    end
    
    NB = ceil(C/Cmax*NBmax);
    
    if NB == 0
        NB = 1;
    end
    
    P1(1,:) = rand(1,4);
    P1(2,:) = round(NB*rand(1,4));

    while sum(P1(2,:)==0)

        P1(1,:) = rand(1,4);
        P1(2,:) = round(NB*rand(1,4));
    
    end   
    
    for i=1:4

        if P1(1,i)<0.5
        
            B1(1,i)=B0(1,i)-P1(2,i)*dds(1,i);
      
        else
     
            B1(1,i)=B0(1,i)+P1(2,i)*dds(1,i);
    
        end
    
        while B1(1,i) > Mxs(1,i) || B1(1,i) < Mns(1,i) 
            
            rnde = rand(1);
            
            if rnde<0.5
        
                B1(1,i)=B0(1,i)-round(NB*rand(1))*dds(1,i);
      
            else
     
                B1(1,i)=B0(1,i)+round(NB*rand(1))*dds(1,i);
    
            end
                
        end 
    
       
        
    end
    
    [Tms,dmp, Gphis, ~] = CoupledParabolicBeam(X,B1(1,3),B1(1,4),B1(1,1),B1(1,2));
    B1(1,5) = OneID(Ga,AcTh,dTh,dt,Tms,dmp,Gphis,0);

    if  B1(1,5)<B0(1,5)
    
        B0 = B1;
        Log = vertcat(Log,B0);
        LogN = LogN+1;
        
    else
        
        BR = B0(1,5)/B1(1,5);
        PT = 1-exp((-1)*C*BR^2);
        
        if rand(1)<PT
           
        B0 = B1;
        B0(1,6) = C;
        Log = vertcat(Log,B0);     
        LogN = LogN+1;
        end 

    end

if Log(LogN,5)<ExitC,break,end    
    
end

end

