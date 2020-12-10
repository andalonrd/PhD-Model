[Gmotion, Path] = uigetfile;
filename = strcat(Path,Gmotion);
A = importdata(filename);

% default plotting colors set to Blue, Red and Green
set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0]);

% Building Parameters

Tn = 0.6;
dmp = 0.05;
alp = 3;
d = 0.4;

Gmotion = (Gmotion);
Tts = strcat( 'T', ' = ', num2str(Tn),'  ');
Tdmp = strcat( '\xi', '=',num2str(dmp), '  ');
Talpha = strcat( '\alpha',' = ' , num2str(alp), '  ');
Tdelta = strcat( '\delta', ' = ', num2str(d), '  ');
LegendId=strcat(Tts,'.',Tdmp,',',Talpha,'.',Tdelta);

[Nts Nx] = size(A);

%Floor Numbers 
Floors = A(1,3:Nx);
% Ordinates of the floors normalized by the building height
Xs = A(2,3:Nx);
% Time ordinates of the acceleration time histories
ts = A(3:Nts,1);
% time interval of the acceleration time histories
dt = ts(4)-ts(3);
% Ground motion First Floor
Gm = A(3:Nts,2);
% Floor Acceleration time histories
FAcTh = A(3:Nts,3:Nx);
Nx = Nx-2;
% Period ordinates of the response spectra (each 0.025s)
Ts =(1:200)*(5/200);
Ts = Ts';
%Response spctra of the Ground motion
GSa = response_spectraA(Gm,0.05,dt,0,0);
% Response to the Building parameters
[IdSa, FSa, IdAcTh ] = LastID(Tn,dmp,alp,d);
% Number of the last floor
RoofN = max(Floors);

% Time histories plotting
figure(1);
subplot(Nx+1,1,1);plot(ts,Gm);
title(strcat( Gmotion ));

for i=1:Nx

    subplot(Nx+1,1,i+1);plot(ts,FAcTh(:,i),ts,IdAcTh(:,i));
    
    if Floors(i) == RoofN
        
        title(strcat('Roof '))
    
    else
        
      title(strcat('Floor ',num2str(Floors(i))));
    
    end

end

legend(LegendId);

figure(2);
subplot(ceil((Nx+1)/2),2,1);plot(Ts,GSa);
title(strcat( Gmotion ));

% Response Spectra Plotting
for i=1:Nx

    subplot(ceil((Nx+1)/2),2,i+1);plot(Ts,FSa(:,i),Ts,IdSa(:,i));
    
    if Floors(i) == RoofN
        
        title(strcat('Roof '))
    
    else
        
      title(strcat('Floor ',num2str(Floors(i))));
    
    end

end

legend('Recorded', LegendId);

