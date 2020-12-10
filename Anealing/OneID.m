function [ IDc ] = OneID(Gm,Tha,Thd,dt,Ts,dmps,Gphis,Wg)

%This function finds the system ID criteria for a single iteration. 

%Input Parameters
%Gm is the ground motion
%St are the Target acceleration Spectra
%dt is the time interval of the samples
%DynP are the dynamic properties of the model being used obtained by the
%function CoupledLegendreBeam
%Wg are the Weights for each time history for the IDcriteria

%Output Parameters
%Id criteria, which is maximized to find the optimal parameters


%Finds the Number of Modes and The number of Locations where Target time
%Histories are available, and the number of samples in the target Th:

[Nx, Nmodes]=size(Gphis);
[ThdR, ~] = size(Thd);
[Wgr, ~] = size(Wg); 
% Defines the default option for the optimization: equal weights for all
% acceleration time histories. 

%if no weights are specified, all floors get the same!
if Wgr == 1
    
    Wa = 0.5;
    Wd = 0.5;
    Wga = 1/Nx*ones(Nx,1);
    Wgd = 1/Nx*ones(Nx,1);
    Wg  = [Wa*Wga Wd*Wgd]; 

end


% Time History Modal Analysis

IDc = 0;

%Note j controls the floor being analyized, i controls the mode being
%combinated. 
for j=1:Nx
   
    Gphi = Gphis(j,:);
    
    Accl = Gm;
    
    if (ThdR ~= 0)
        
        disp = 0*Gm;
  
    end
    
    for i=1:Nmodes
        
        iAccl = Gphi(i)*sdofrhaA(Gm,Ts(i), dmps(i),dt,0,0);
        
        % only if the mode constributes to a 2% of the maximum: (except
        % the first one that altways goes!
        
        if ( i==1 || max(abs(iAccl))>0.02*max(abs(Accl)))
            
            Accl = Accl + iAccl;
        
        end
        
        % If Displacement time histories are being considered:
        if ( ThdR ~= 0)
            
            idisp = Gphi(i)*sdofrhaD(Gm,Ts(i),dmps(i),dt,0,0);
            
            if (i==1 || max(abs(idisp))>0.02*max(abs(disp)))
                
                 disp = disp + idisp;
            
            end
            
        end
        
    end
    
    % If Displacement time histories are being considered:
    if (ThdR ~= 0)
        
        EF = max(abs(Tha(:,j)))/max(abs(Thd(:,j)));
        IDc = IDc + Wg(j,1)*sum((Accl-Tha(:,j)).^2)+Wg(j,2)*EF*sum((disp-Thd(:,j)).^2);
    
    else
        
        IDc = IDc +2* Wg(j,1)*sum((Accl-Tha(:,j)).^2);
    
    end    
end


end

