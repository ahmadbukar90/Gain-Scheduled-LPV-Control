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
% G1 = ss(A,B(:,1), C, D(:,1));
Gd = tf(1000/Iz);
%% Define the time-varying real parameter.
Ur = 1; % rough_road
d_max = 0.40;%degrees
del = pgrid('del',linspace(-d_max,d_max,11));
%rho = pgrid('rho',-d_max:0.5:d_max)
del.RateBounds = [-0.6 0.6];

v_min = 10; %m/s
v_max = 40;
v = pgrid('v',linspace(v_min,v_max,11));
v.RateBounds = [5 50];
% Define the A, B, C, and D matrices of the LPV system in Equation (1)
%%
%v=10;
p1 = 1/v;
p2 = 1/(v.*v);
p3 = cos(del);
% p3=1;

%%
A11 = p1*Ur*(-C_f-C_r)/m;
A12 = 1+ p2*Ur*(l_r*C_r-l_f*C_f)/m;
A21 = Ur*(l_r*C_r-l_f*C_f*p3)/Iz;
A22 = p1*Ur*(-l_f*l_f*C_f*p3-l_r*l_r*C_r)/Iz;

B11 = p1* C_f/m;
B21 = p3* l_f*C_f/Iz;
B22 = S_rr*R*t_r/2/Iz;
B23 = -B22;

A = [A11 A12;A21 A22];
B = [B11 0 0;B21 B22 B23];
C = [1 0];
D=[0 0 0];

% Form the grid-based parameter-varying state-space system:
G = pss(A,B,C,D);

%% open loop simulation
t = linspace(0,10,101);
ptraj.time = t;%linspace(0,20,500);
ptraj.v = repmat(10,[101 1]);
ptraj.del = 0.3*sin(t);

figure (1)
plot(t,ptraj.v)
hold on
plot(t,ptraj.del)

figure (2)
lpvstep(G,ptraj); %LPV simulation response

H_inf_norm = lpvmax( norm(G,inf));% Find peak Hinf over domain
Lpv_norm = lpvnorm(G);%Bound induced L2 norm over allowable (rate unbounded) trajectories

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
W3 = Wtbrj;
W4 = Wtbrj;
% bode(W1)
We = W1;
Wu = [W2 0 0;0 W3 0;0 0 W4];
% Wind effect
Gd = tf(1000/Iz);

%% Generalized plant P is found with function sysic

systemnames='G Gd We Wu';
input_to_G='[u]';
input_to_Gd='[d]';
input_to_We='[r-G-Gd]';
input_to_Wu='[u]';
inputvar='[r;d;u(3)]';
outputvar='[We;Wu;r-G-Gd]';
sysoutname   = 'P';
cleanupsysic = 'yes';

sysic;

%%
% Basis function,
b1 = basis(1,0);
b2 = basis(p1,'v',1);
b3 = basis(p2,'v',1);
b4 = basis(p3,'del',1);

% basis(delta,'delta',1);
% bp1 = basis(pcos,'rho',-psin);
% bsin = basis(psin,'rho',pcos);
Xb = [b1,b2,b3,b4];
Yb = Xb;

nmeas=1;
ncon=3;
% %% grid Hinf control design
% [Khinf,normhinf,info] = hinfsyn(P,nmeas,ncon);
% peakhinf = norm(info.Data(:),inf);

%% grid LPV non rate-Bounded Control Design
opt = lpvsynOptions('BackOffFactor',1.02);
[Klpvn,normlpvn,info] = lpvsyn(P,nmeas,ncon,opt);
%% grid LPV Rate-Bounded Control Design
opt = lpvsynOptions('BackOffFactor',1.02);
[Klpv,normlpv,info] = lpvsyn(P,nmeas,ncon,Xb,Yb,opt);

% Eliminate rate dependence of LPV Rate-Bounded Controller
klpv2r = lpvinterp(Klpv,{'delDot'},{0});
Klpve = lpvelimiv(klpv2r);
%% norm comparison

% normcomp = [peakhinf normlpvn normlpv];

%% Basis function derivative,
b1 = basis(1,0);
b2 = basis(p1,'v',1);
b3 = basis(p3,'del',1);

% basis(delta,'delta',1);
% bp1 = basis(pcos,'rho',-psin);
% bsin = basis(psin,'rho',pcos);
Xb = [b1,b2,b2^2,b3,b3^2];
Yb = Xb;

nmeas=1;
ncon=3;
% rb = logspace(-2,log10(4),15);
% for i=1:numel(rb)
%   P.Parameter.del.RateBounds = [-rb(i) rb(i)];
%   bnds4(i) = lpvnorm(P,Xb(1:4));
%   bnds1(i) = lpvnorm(P,Xb(1:2));
%   bnds2(i) = lpvnorm(P,Xb(1:3));
%   bnds3(i) = lpvnorm(P,Xb);
% end

%% grid LPV Rate-Bounded Control Design
opt = lpvsynOptions('BackOffFactor',1.02);
[Klpvd,normlpv,info] = lpvsyn(P,nmeas,ncon,Xb,Yb,opt);

% Eliminate rate dependence of LPV Rate-Bounded Controller
klpv2r = lpvinterp(Klpv,{'delDot'},{0});
Klpvde = lpvelimiv(klpv2r);
% % Generate LPV controller for different rate bounds
% rb = [0.1 1:10]; 
% for ii = 1:numel(rb);
%      P.Parameter.del.RateBounds = [-rb(ii) rb(ii)];
%      [~,nl1] = lpvsyn(P,nmeas,ncon,Xb,Yb,opt);
%      nlpv3(ii) = nl1;
% end

%%
omega=logspace(-5,5,10000);

%clrb_grid = lft(P,Klpv);  % Rated bound closed loop system
clrbe_grid = lft(P,Klpve);  % eleiminating the rated bound
clnrn_grid = lft(P,Klpvn); % Non rated bound closed loop system

figure (3)
subplot(2,1,1);
sigma(clrbe_grid(1,1)/W1,'r',clnrn_grid(1,1)/W1,'b',1/W1,'k--',omega);
title('Ref Sensitivity: Z1/r')
legend('S-RTB','S-NRTB','1/We')

subplot(2,1,2);
sigma(clrbe_grid(1,2)/W1,'r',clnrn_grid(1,2)/W1,'b',1/W1,'k--',omega);
title('Wind effect sensitivity: Z1/d')
legend('SG-RTB','SG-NRTB','1/We')
%%
figure (4)
subplot(2,1,1);
sigma(clrbe_grid(2,1)/W2,'r',clnrn_grid(2,1)/W2,'b',1/W2,'k--',omega);
title('Stearing control sensitivity: Z2/r')
legend('KS-RTB','KS-NRTB','1/Wu1')

subplot(2,1,2);
sigma(clrbe_grid(2,2)/W2,'r',clnrn_grid(2,2)/W2,'b',1/W2,'k--',omega);
title('Stearing control sensitivity: Z2/d')
legend('KS-RTB','KS-NRTB','1/Wu1')

figure (5)
subplot(2,1,1);
sigma(clrbe_grid(3,1)/W2,'r',clnrn_grid(3,1)/W3,'b',1/W3,'k--',omega);
title('Braking control sensitivity: Z3/r')
legend('KS-RTB','KS-NRTB''1/Wu1')

subplot(2,1,2);
sigma(clrbe_grid(3,2)/W2,'r',clnrn_grid(3,2)/W3,'b',1/W3,'k--',omega);
title('Braking control sensitivity: Z3/d')
legend('KS-brk','1/Wu2')
%%
omega=logspace(-5,5,10000);

% clrb_grid = lft(P,Klpvd);  % Rated bound closed loop system
clrbe_grid = lft(P,Klpvde);  % eleiminating the rated bound
clnrn_grid = lft(P,Klpvn); % Non rated bound closed loop system

figure (6)
subplot(2,1,1);
sigma(clrbe_grid(1,1)/W1,'r',clnrn_grid(1,1)/W1,'b',1/W1,'k--',omega);
title('Ref Sensitivity: Z1/r')
legend('S-RTB','S-NRTB','1/We')

subplot(2,1,2);
sigma(clrbe_grid(1,2)/W1,'r',clnrn_grid(1,2)/W1,'b',1/W1,'k--',omega);
title('Wind effect sensitivity: Z1/d')
legend('SG-RTB','SG-NRTB','1/We')

figure (7)
subplot(2,1,1);
sigma(clrbe_grid(2,1)/W2,'r',clnrn_grid(2,1)/W2,'b',1/W2,'k--',omega);
title('Stearing control sensitivity: Z2/r')
legend('KS-RTB','KS-NRTB','1/Wu1')

subplot(2,1,2);
sigma(clrbe_grid(2,2)/W2,'r',clnrn_grid(2,2)/W2,'b',1/W2,'k--',omega);
title('Stearing control sensitivity: Z2/d')
legend('KS-RTB','KS-NRTB','1/Wu1')

figure (8)
subplot(2,1,1);
sigma(clrbe_grid(3,1)/W2,'r',clnrn_grid(3,1)/W3,'b',1/W3,'k--',omega);
title('Braking control sensitivity: Z3/r')
legend('KS-RTB','KS-NRTB''1/Wu2')

subplot(2,1,2);
sigma(clrbe_grid(3,2)/W2,'r',clnrn_grid(3,2)/W3,'b',1/W3,'k--',omega);
title('Braking control sensitivity: Z3/d')
legend('KS-RTB','KS-NRTB''1/Wu2')