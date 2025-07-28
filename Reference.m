function xDot = Reference(U)

parameters

Ur = 0.5;
beta = U(1);
del = U(2);
v = U(3);
psid = U(4);

betadot = (Ur*(-C_f-C_r)/(m*v))*beta + (1+ Ur*(l_r*C_r-l_f*C_f)/(m*v^2))*psid + (C_f/(m*v))*del;
psiddot = (Ur*(l_r*C_r-l_f*C_f)/Iz)*beta + (Ur*(-l_f*l_f*C_f-l_r*l_r*C_r)/(Iz*v))*psid + (l_f*C_f/Iz)*del; 

xDot = [betadot psiddot];
end