function [ A ] = CoupledUniformBeam(X,alp,Tn,dmp)

B = importdata('Gvalues.txt');

[~, Nmodes] = size(B);
dalp = B(2,1)-B(1,1);
Nv = floor(alp/dalp)+1;
gmv = B(Nv,2:Nmodes)';
clear('B');
Nmodes =Nmodes-1;
Nx = length(X);

if X(Nx)==1

    X(Nx)=0.999999;
    
end

A = zeros(2+Nx,5);


if length(dmp) == 1

    dmp = ones(5,1)*dmp;

end

for j=1:Nmodes
   
    gm=gmv(j);
    bt = sqrt(alp^2+gm^2);
    nm = (bt*gm^2*sin(gm)-gm^3*cos(gm)-gm*bt^2*exp(-bt))/(bt*gm^2*cos(gm)+bt^3*cosh(bt))+gm/bt;
    f1 = exp(bt)*(nm-gm/bt);
    f2 = exp(-bt)*(nm+gm/bt);
    MPFN = 1/(2*bt^2*gm)*(2*bt^2+2*gm^2-2*bt^2*cos(gm)-2*bt^2*nm*sin(gm)+bt*gm*(f1-f2));
    MPFD = 0.5-1/(2*gm)*sin(gm)*cos(gm)-nm/gm*sin(gm)^2+nm^2/2+nm^2/(2*gm)*sin(gm)*cos(gm)+bt/(bt^2+gm^2)*(f1-f2)*(sin(gm)-nm*cos(gm))-gm/((gm^2+bt^2))*(f1+f2)*(cos(gm)+nm*sin(gm))+sinh(bt)/(4*bt)*(f1*(nm-gm/bt)+f2*(nm+gm/bt))+f1*f2/2;
    MPF = MPFN/MPFD; 
    A(1,j)=Tn*gmv(1)/gm*sqrt((gmv(1)^2+alp^2)/(gm^2+alp^2));
    A(2,j)=dmp(j);
    
    for k=1:Nx
        
        x=X(k);        
        A(2+k,j) = MPF*(sin(gm*x)-nm*cos(gm*x)+0.5*(nm-gm/bt)*exp(bt*x)+0.5*(nm+gm/bt)*exp(-bt*x));
    
    end
    

end

