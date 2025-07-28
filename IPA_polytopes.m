% IPA Project, ALI and Ahmad  
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
%  x = [beta psidt];
%  u = [delta_f Tbl Tbr];
%  y = [beta psiddt];
parameters
%% Reference model = psidt_ref
Pr = referencemodel(m,Iz,l_f,l_r,t_f,t_r,C_f,C_r,S_r,v);
%% Define the time-varying real parameter.
% delta=[-30,30]degrees, Vel=[6,50]

meo = 1;
d_min = 0;%rad
d_max = 0.6;%rad
v_min = 6; %m/s
v_max = 40;%m/s

% normalize values of P1=1/v
p1min = 1/v_max;
p1max = 1/v_min;

% normalize values of P2=1/(v*v)
p2min = 1/(v_max.*v_max);
p2max = 1/(v_min.*v_min);

% normalize values of P3 = del 
p3max = cos(d_min);
p3min = cos(d_max);

%% ===== low pass Filter ===========
wf=40;
Af= -wf;
Bf= wf;
Cf= 1;
Df= 0;
F=ss(Af,Bf,Cf,Df);
%% sub- blocks

A02 = Af;

A10 =  C_f/m; % B11 = p1* C_f/m;
A11 = meo*(-C_f-C_r)/m;
A12 = meo*(l_r*C_r-l_f*C_f)/m;

A20 =  l_f*C_f/Iz; % B21 = l_f*C_f/Iz;
% A21 = (l_r*C_r-l_f*C_f*p3)/Iz;
% A22 = (-l_f*l_f*C_f*p3-l_r*l_r*C_r)/Iz;

%%
A000 = [A02 0 0;p1min*A10 p1min*A11 (1+ p2min*A12);
        p3min*A20 meo*(l_r*C_r-l_f*C_f*p3min)/Iz p1min*meo*(-l_f*l_f*C_f*p3min-l_r*l_r*C_r)/Iz]; % p1min,p2min,p3min

A001 = [A02 0 0;p1min*A10 p1min*A11 (1+ p2min*A12);
        p3max*A20 meo*(l_r*C_r-l_f*C_f*p3max)/Iz p1min*meo*(-l_f*l_f*C_f*p3max-l_r*l_r*C_r)/Iz]; % p1min,p2min,p3max
   
A010 = [A02 0 0;p1min*A10 p1min*A11 (1+ p2max*A12);
        p3min*A20 meo*(l_r*C_r-l_f*C_f*p3min)/Iz p1min*meo*(-l_f*l_f*C_f*p3min-l_r*l_r*C_r)/Iz]; % p1min,p2max,p3min
    
A011 = [A02 0 0;p1min*A10 p1min*A11 (1+ p2max*A12);
        p3max*A20 meo*(l_r*C_r-l_f*C_f*p3max)/Iz p1min*meo*(-l_f*l_f*C_f*p3max-l_r*l_r*C_r)/Iz]; % p1min,p2max,p3max
    
A100 = [A02 0 0;p1max*A10 p1max*A11 (1+ p2min*A12);
        p3min*A20 meo*(l_r*C_r-l_f*C_f*p3min)/Iz p1max*meo*(-l_f*l_f*C_f*p3min-l_r*l_r*C_r)/Iz]; % p1max,p2min,p3min
    
A101 = [A02 0 0;p1max*A10 p1max*A11 (1+ p2min*A12);
        p3max*A20 meo*(l_r*C_r-l_f*C_f*p3max)/Iz p1max*meo*(-l_f*l_f*C_f*p3max-l_r*l_r*C_r)/Iz]; % p1max,p2min,p3max
    
A110 = [A02 0 0;p1max*A10 p1max*A11 (1+ p2max*A12);
        p3min*A20 meo*(l_r*C_r-l_f*C_f*p3min)/Iz p1max*meo*(-l_f*l_f*C_f*p3min-l_r*l_r*C_r)/Iz]; % p1max,p2max,p3min
    
A111 = [A02 0 0;p1max*A10 p1max*A11 (1+ p2max*A12);
        p3max*A20 meo*(l_r*C_r-l_f*C_f*p3max)/Iz p1max*meo*(-l_f*l_f*C_f*p3max-l_r*l_r*C_r)/Iz]; % p1max,p2min,p3max

%% testing stability at vertex

%%% Create variables matrix
X  = sdpvar(3,3,'symmetric'); 

epsi = eps;
%%% LMIs definition 
H = [X ];

M1 = A000'*X+X*A000;
H1  = [M1];

M2 = A001'*X+X*A001;
H2  = [M2];

M3 = A010'*X+X*A010;
H3  = [M3];

M4 = A011'*X+X*A011;
H4  = [M4];

M5 = A100'*X+X*A100;
H5  = [M5];

M6 = A101'*X+X*A101;
H6  = [M6];

M7 = A110'*X+X*A110;
H7  = [M7];

M8 = A111'*X+X*A111;
H8  = [M8];

F  = [H1<0,H2<0,H3<0,H4<0,H5<0,H6<0,H7<0,H8<0,H>0];

%%% Find feasible solution 
%ops      = sdpsettings('solver',sedumi);
solution = solvesdp(F,[]);

%solution = solvesdp(F);

X  = double(X)
%%% Check solution
VPx = eig(X)
for i=1:size(X,1)
    if (VPx(i) <= 0)
        disp('Error in the Lyapunov function')
    end;
end;
checkset(F)

%% Stability and Open loop analysis
 
figure (1)  %stability analysis
plot(VPx,'* r')
title('Eigenvalues of lyaponov function')
xlabel('Index')
ylabel('Eigenvalue')
%%
B00 = Bf;
B22 = S_rr*R*t_r/2/Iz;
B23 = -B22;

B = [B00 0 0;0 0 0;0 B22 B23];
C = [0 0 1];
D = [0 0 0];

G000=ss(A000,B,C,D);
G001=ss(A001,B,C,D);
G010=ss(A010,B,C,D);
G011=ss(A011,B,C,D);
G100=ss(A100,B,C,D);
G101=ss(A101,B,C,D);
G110=ss(A110,B,C,D);
G111=ss(A111,B,C,D);

%%
t = 0:0.01:5;
Tb =zeros(length(t),1);
del =zeros(length(t),1);

for x =1:1:length(t)
    del(x)=10;
end
v = [del Tb Tb];
y = lsim(G000,v,t);
plot(t,y)
%% Yaw Parameters
Slope = 0.2;
angle = 40;
Tstart = 0.1;
double = 1.0;
%% Weightings for Robust Controller

% Performances
Ge      = .1;
w1      = 2*pi*1;
w2      = w1/Ge;
WeYawdt = tf([1/w2 1],[1/w1 1])/(2*Ge); 

% Control
alpha   = 100;
w5      = 2*pi*1;%1
w6      = 2*pi*10;%10
delta   = (w6+w5)/2;
Wdelta  = 1*tf(conv([1/w5 1],[1/w6 1]),conv([alpha*w6 1],[1/(alpha*w6) 1]));
Gdelta  = bode(Wdelta,delta);
Wdelta  = 5e-3*inv(Gdelta)*Wdelta;%-2
w7      = 2*pi*10;
Wtbrj    = 1e-4*tf([1/w7 1],[1/(alpha*w7) 1]);%-4
WtbL   = Wtbrj;
WtbR   = Wtbrj;

W1 = WeYawdt;
W2 = Wdelta;
W3 = WtbL;
W4 = WtbR;
%
WactS   = 1*tf(1,[1/w6 1]);
WactBl  = 1*tf(1,[1/w7 1]);
WactBr  = WactBl;
We = W1;
Wu = [W2 0 0;0 W3 0;0 0 W4];
% Wind effect
Gd = tf(1000/Iz);
%% Generalized plant P is found with function sysic
G = G000;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P000';
cleanupsysic = 'yes';

sysic;

G = G001;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P001';
cleanupsysic = 'yes';

sysic;

G = G010;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P010';
cleanupsysic = 'yes';

sysic;

G = G011;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P011';
cleanupsysic = 'yes';

sysic;

G = G100;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P100';
cleanupsysic = 'yes';

sysic;

G = G101;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P101';
cleanupsysic = 'yes';

sysic;

G = G110;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P110';
cleanupsysic = 'yes';

sysic;

G = G111;

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P111';
cleanupsysic = 'yes';

sysic;

listP={P001,P001,P010,P011,P100,P101,P110,P111};
nmeas=1;
ncon=3;
sp=1;
percentage=0;

%% LPV/hinf controller synthesis - One Degree of freedom
 
[listK,listCL,gopt] = lmiHinfPolytope(listP,nmeas,ncon,sp,percentage,'sedumi')
omega=logspace(-5,5,10000);

%% sensitivity functions
% K1
figure(2)

subplot(2,1,1);
sigma(listCL{1}(1,1)/W1,listCL{2}(1,1)/W1,listCL{3}(1,1)/W1,listCL{4}(1,1)/W1,listCL{5}(1,1)/W1,listCL{6}(1,1)/W1,listCL{7}(1,1)/W1,listCL{8}(1,1)/W1,1/W1,'k--',omega)
title('Ref Sensitivity: Z1/r')
legend('S','1/We')

subplot(2,1,2);
sigma(listCL{1}(1,2)/W1,listCL{2}(1,2)/W1,listCL{3}(1,2)/W1,listCL{4}(1,2)/W1,listCL{5}(1,2)/W1,listCL{6}(1,2)/W1,listCL{7}(1,2)/W1,listCL{8}(1,2)/W1,1/W1,'k--',omega)
title('Wind effect sensitivity: Z1/d')
legend('SG','1/We')

figure (3)

subplot(2,1,1);
sigma(listCL{1}(2,1)/W2,listCL{2}(2,1)/W2,listCL{3}(2,1)/W2,listCL{4}(2,1)/W2,listCL{5}(2,1)/W2,listCL{6}(2,1)/W2,listCL{7}(2,1)/W2,listCL{8}(2,1)/W2,1/W2,'k--',omega)
title('Stearing control sensitivity: Z2/r')
legend('KS-del','1/Wu1')
subplot(2,1,2);
sigma(listCL{1}(2,2)/W2,listCL{2}(2,2)/W2,listCL{3}(2,2)/W2,listCL{4}(2,2)/W2,listCL{5}(2,2)/W2,listCL{6}(2,2)/W2,listCL{7}(2,2)/W2,listCL{8}(2,2)/W2,1/W2,'k--',omega)
title('Stearing control sensitivity: Z2/d')
legend('KS-del','1/Wu1')

figure (4)

subplot(2,1,1);
sigma(listCL{1}(3,1)/W3,listCL{2}(3,1)/W3,listCL{3}(3,1)/W3,listCL{4}(3,1)/W3,listCL{5}(3,1)/W3,listCL{6}(3,1)/W3,listCL{7}(3,1)/W3,listCL{8}(3,1)/W3,1/W3,'k--',omega)
title('Braking control sensitivity: Z3/r')
legend('KS-brk','1/Wu2')
subplot(2,1,2);
sigma(listCL{1}(3,2)/W3,listCL{2}(3,2)/W3,listCL{3}(3,2)/W3,listCL{4}(3,2)/W3,listCL{5}(3,2)/W3,listCL{6}(3,2)/W3,listCL{7}(3,2)/W3,listCL{8}(3,2)/W3,1/W3,'k--',omega)
title('Braking control sensitivity: Z3/d')
legend('KS-brk','1/Wu2')

