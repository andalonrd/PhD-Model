function [FiD ] = ID(T0,dmp0,alp0,d0 )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[BaseMotion, Path] = uigetfile;
warning('off','MATLAB:nearlySingularMatrix');
filename = strcat(Path,BaseMotion);

InFile = importdata(filename);
[Nth Nx] = size(InFile);
X = InFile(2,3:Nx);
Ga = InFile(3:Nth,2);
AcTh = InFile(3:Nth,3:Nx);
dt = InFile(4,1)-InFile(3,1);
clear ('InFile');
FiD=OneID(Ga,AcTh,0,dt,CoupledLegendreBeam(X,alp0,d0,T0,dmp0),0);
FiD = 1/FiD;
end

