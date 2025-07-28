function xDot = Nonbicycle(U)

parameters

Ur = 1;
beta = U(1);
del = U(2);
Tbrf = U(3);
Tbrr = U(4);
Mdz = U(5);
v = U(6);
psid = U(7);

betadot = (Ur*(-C_f-C_r)/(m*v))*beta + (1+ Ur*(l_r*C_r-l_f*C_f)/(m*v^2))*psid + (C_f/(m*v))*del;
psiddot = (Ur*(l_r*C_r-l_f*C_f*cos(del))/Iz)*beta + (Ur*(-l_f*l_f*C_f*cos(del)-l_r*l_r*C_r)/(Iz*v))*psid + (l_f*C_f*cos(del)/Iz)*del +(1/Iz)*Mdz + (S_rr*R*t_r/2/Iz)*Tbrf - (S_rr*R*t_r/2/Iz)*Tbrr; 

xDot = [betadot psiddot];
end