function u   = MPC_control(x1,x2,xref)

parameters

Rad   = .3;
v= 26;    % Lateral Velocity

N =20; % Prediction horizon
ud =[ 0;0;0 ]; % Desired input
Q = diag([1 10]);
R = diag([1 1 1]);
P=10*eye(2); %weight matrix for terminal states
nu=3;
%%%%-------Tracking Trajectory------%%%%    
xd = [0 ; xref];
x=[x1;x2];
x0 = x;


A_c = [(-C_f-C_r)/(m*v) 1+(l_r*C_r-l_f*C_f)/(m*v*v);
     (l_r*C_r-l_f*C_f)/Iz (-l_f*l_f*C_f-l_r*l_r*C_r)/(v*Iz)];
B_c = [C_f/(m*v) 0 0;
     l_f*C_f/Iz S_rr*Rad*t_r/2/Iz -S_rr*Rad*t_r/2/Iz];
%C = [0 1];
%D=[0 0 0];
% Note that we are assuming a full state feedback
%%%%%%%% --------------------------%%%%%%%%%%%%%%%




%%%%%%--------Convert continous time dynamics to discrete time dynamics ---%%%%%  
tau = 0.2; %sampling time
[A_d,B_d] = c2d(A_c,B_c,tau); %discretization

x_max = [inf; inf];   x_min = [-inf; -inf];
%x_max = [ inf; inf];   x_min = [-inf; -inf];
% Constraints on the inputs
 u_max = [inf; inf; inf ];       u_min = [-inf;-inf;-inf];
%%%%-----Compute the cost------%%%%%%
[H,h,c] = compute_cost_matrices(Q,R,P,N,xd,ud);

%%%%-----Equality constraints------%%%%%%
[Phi,Psi] =construct_Gamma(A_d,B_d,N);
A_eq = Psi;
B_eq= -1*(Phi*x0);

%%%%-----In-equality constraints------%%%%%%
[A_ineq, B_ineq] = constuct_constraint_mat(u_max, u_min, x_max, x_min, N);

%%%%%------Solving the optimization problem ------%%%%%
solution = quadprog(2*H,h,A_ineq,B_ineq,A_eq,B_eq);
u=solution(1:nu)'; %computing the optimal input
%lesx(i+1,:)=A_d*(lesx(i,:))'+B_d*lesu(i,:)'; % state update rule

u=u;
