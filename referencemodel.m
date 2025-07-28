% Author: Charles Poussot-Vassal
% 
% Description
% Function that compute the bicycle vehicle linear model as a function of
% the speed. This function is usefull for controller synthesis at a defined
% speed. Note that steer angle input is in degree and speed in m/s. Braking
% here only affects the yaw dynamic.
% 
% Input
%  m   : vehicle total mass
%  Ix  : vehicle roll inertia 
%  Iz  : vehicle yaw inertia 
%  l_f : distance COG front
%  l_r : distance COG rear
%  h   : height of the vehicle COG
%  C_f : linearized lateral friction coeff of front tire
%  C_f : linearized lateral friction coeff of rear tire
%  v   : longitudinal speed of the vehicle [m/s]
% 
% Output
%  P : system state space
% 
%  x = [beta psidt]
%  u = [delta_f Mdz Tbl Tbr]
%  y = [betadt beta psiddt psidt]
%
% bicycle = bicycleLinearModel0(m,Iz,l_f,l_r,t_f,t_r,C_f,C_r,S_r,v)

function Pr = referencemodel(m,Iz,l_f,l_r,t_f,t_r,C_f,C_r,S_r,v)
parameters

% v = 25;m/s

Ar = [(-C_f-C_r)/(m*v) 1+(l_r*C_r-l_f*C_f)/(m*v*v);
     (l_r*C_r-l_f*C_f)/Iz (-l_f*l_f*C_f-l_r*l_r*C_r)/(v*Iz)];
Br = [C_f/(m*v);l_f*C_f/Iz ];
Cr = [0 1];
Dr = [0];
Pr = ss(Ar,Br,Cr,Dr);