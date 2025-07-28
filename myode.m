
function xDot = myode(t,x,U,delt,del)

xDot = zeros(2,1);
parameters

% del = interp1(delt,del,t);

Ur = 0.5;
beta = x(1);
% del = U(1);
Tbrf = U(1);
Tbrr = U(2);
v = U(3);
psid = x(2);

xDot(1) = (Ur*(-C_f-C_r)/m*v)*beta + (1+ Ur*(l_r*C_r-l_f*C_f)/m*v^2)*psid + (C_f/m)*del;
xDot(2) = (Ur*(l_r*C_r-l_f*C_f*cos(del))/Iz)*beta + (Ur*(-l_f*l_f*C_f*cos(del)-l_r*l_r*C_r)/Iz)*psid + (l_f*C_f*cos(del)/Iz)*del + (S_rr*R*t_r/2/Iz)*Tbrf - (S_rr*R*t_r/2/Iz)*Tbrr; 
end