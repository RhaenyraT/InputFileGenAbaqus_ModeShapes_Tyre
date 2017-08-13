function [qmem_P0pos]  =NL_sidewallpos(sinput)

L          = zeros(length(sinput.urpos),1);
w          = zeros(length(sinput.urpos),1);
qmem_P0pos = zeros(length(sinput.urpos),1);
w(1)       = sinput.w; %0.8994; 

for i= 2:length(sinput.urpos)
    
    L(i)          =   sinput.Tire_Belt_Radius - sinput.Tire_Rim_Radius + sinput.urpos(i);
    
    w(i)          =  fsolve(@(w)(sin(w/2)-L(i)*w/(2*sinput.lstring)),w(i-1));
    
    qmem_P0pos(i) =   -(sinput.lstring/w(i))*(cos(w(i)/2));

end
