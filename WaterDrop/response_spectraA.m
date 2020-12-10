% Linear elastic response response spectra
%
% Written by Eduardo Miranda and marios Kyriakides on 10/27/05
%

function [sa]=response_spectraA(accel,xi,dt,d0,v0)

% initializes vectors where response parameters will be stored
t1=zeros(200,1);
sa=zeros(200,1);
sv=zeros(200,1);
sd=zeros(200,1);
psa=zeros(200,1);
psv=zeros(200,1);

for i=1:200;
    t1(i)=i*0.025;
    % Performs response history analysis using recursive method
    [rd,rv,ra,aa,f]=sdofrha(accel,t1(i),xi,dt,d0,v0);

    % Computes peaks
    if max(max(aa),max(abs(aa)))>sa(i); 
        sa(i)=max(max(aa),max(abs(aa)));
    end
    if max(max(rv),max(abs(rv)))>sv(i);
        sv(i)=max(max(rv),max(abs(rv)));
    end
    if max(max(rd),max(abs(rd)))>sd(i);
        sd(i)=max(max(rd),max(abs(rd)));
    end
    psa(i)=4*pi*pi*sd(i)/t1(i)/t1(i);
    psv(i)=2*pi*sd(i)/t1(i);
end

end
