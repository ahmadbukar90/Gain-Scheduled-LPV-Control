function [A_ineq, B_ineq] = constuct_constraint_mat(u_max, u_min, x_max, x_min, N)

nx = length(x_min);
nu = length(u_max);

rows_Gk = 2*nu+2*nx; % rows of G_k matrix
cols_Gk = nx+nu; % columns of G_k matrix

G_k = zeros(rows_Gk,cols_Gk);
d_k = zeros(rows_Gk,1);

%%%-----Constructing G_k --------%%
G_k(1:nu,1:nu) = eye(nu);
G_k(nu+1:2*nu,1:nu) = -1*eye(nu);
G_k(2*nu+1:(2*nu)+nx,nu+1:end) = eye(nx);
G_k((2*nu)+nx+1:end,nu+1:end) = -1*eye(nx);

%%%-----Constructing d_k --------%%
d_k = [u_max;-1*u_min;x_max;-1*x_min];

%%%-----Constructing A_inq and B_inq --------%%
A_ineq =[];
B_ineq= [];
temp_Vec = zeros(1,N); %temporary variable for constructing S_k

for k =1:N
   
   %%%%-----Computing S_k----------%%%
   temp_var = temp_Vec;
   temp_var(1,k)=1;
   S_k = otimes(temp_var,eye(nx+nu));% Computing Sk    
   %%%%%%%----------------%%%%%%%%
   
   A_ineq = [A_ineq ; G_k*S_k];
   B_ineq = [B_ineq ; d_k];
   
end

end