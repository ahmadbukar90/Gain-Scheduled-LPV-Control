function xDot = odefcn(t,x,U)

parameters

Ur = 0.5;
beta = x(1);
del = U(1);
Tbrf = U(2);
Tbrr = U(3);
Mdz = U(4);
v = U(5);
psid = x(2);

betadot = (Ur*(-C_f-C_r)/m*v)*beta + (1+ Ur*(l_r*C_r-l_f*C_f)/m*v^2)*psid + (C_f/m)*del;
psiddot = (Ur*(l_r*C_r-l_f*C_f*cos(del))/Iz)*beta + (Ur*(-l_f*l_f*C_f*cos(del)-l_r*l_r*C_r)/Iz)*psid + (l_f*C_f*cos(del)/Iz)*del +(1/Iz)*Mdz + (S_rr*R*t_r/2/Iz)*Tbrf - (S_rr*R*t_r/2/Iz)*Tbrr; 

xDot = [betadot;psiddot];
end