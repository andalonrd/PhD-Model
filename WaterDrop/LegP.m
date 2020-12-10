

function [ h ] = LegP(v,x )

v0 = v-floor(v);
v1 = v0 + 1;
h0 = LegPb(v0,x);
h1 = LegPb(v1,x);
h2 = ((2*v0+3)*x*h1-(v0+1)*h0)/(v0+2);
vi = 2 +v0;

if v< 2

    h = LegPb(v,x);
    
    else
        
    while vi ~= v
   
        v0 = v0+1;
        h0 = h1;
        h1 = h2;
        h2 = ((2*v0+3)*x*h1-(v0+1)*h0)/(v0+2);
        vi = v0+2;
    end
    
    h=h2;

end


end

