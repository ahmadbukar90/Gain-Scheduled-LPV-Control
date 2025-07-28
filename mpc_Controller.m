function [mpc_control] =mpc_Controller(A,B,N,Q,R,P,xd,ud,x)

[nx,nu] = size(B);
x =x';

%%%%%-----Constructing Gamma and Cost Matrices-------%%%%%%%

[Phi,Psi] =construct_Gamma(A,B,N);
[H,h,c] = compute_cost_matrices(Q,R,P,N,xd,ud);

%%%%------Computing the optimal input-------%%%%%%%

temp0 = -1*[h ;Phi*x]; %% -[h;Phix ]
temp1 = [2*H Psi']; %% [2H Psi']
temp2 = [Psi zeros(size(Psi,1),size(temp1,2)-size(Psi,2))]; %% [Psi 0]
temp = vertcat(temp1,temp2); %%[2H Psi'; Psi 0]

solution = inv(temp)*temp0; %%[z* lamda*]
mpc_control = solution(1:nu); %uo*

end