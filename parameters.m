%% System parameters
l_f    = 1.4;
l_r    = 1;
t_f    = 1.4;
t_r    = 1.4;
h      = .4;
ms     = 350*4;
mus_fl = 35;
mus_fr = 35;
mus_rl = 32.5;
mus_rr = 32.5;
mus_f  = mus_fl+mus_fr;
mus_r  = mus_rl+mus_rr;
m      = ms+mus_f+mus_r;
m_f    = ms*l_f/(l_f+l_r) + mus_f;
m_r    = ms*l_r/(l_f+l_r) + mus_r;
R      = .3;
Iw     = 1;
g      = 9.81;
Ix     = 250;
Iy     = 1400;
Iz     = m*l_r*l_f;


C_f    = 84085/2; %63000*mus_f/2;
C_r    = 87342/2; %63000*mus_r/2;%20000;%

m_f    = ms*l_f/(l_f+l_r) + mus_f;
m_r    = ms*l_r/(l_f+l_r) + mus_r;
S_f    = m_f*g;
S_r    = m_r*g;

%P = bicycleLinearModel0(m,Iz,l_f,l_r,t_f,t_r,C_f,C_r,S_r,v)
R   = .3;
S_rr = S_r/R;

 v = 15 ;% m/s
