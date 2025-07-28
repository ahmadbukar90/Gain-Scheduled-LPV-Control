function [H,h,c] = compute_cost_matrices(Q,R,P,N,xd,ud)

temp_Vec = zeros(1,N); % Temporary variable that helps in computing Sk
nx = size(Q,1); %system matrix rows
nu = size(R,2); %input matrix rows

for(k =1:N-1)
    
temp_var = temp_Vec;
temp_var(1,k)=1;
S_k = otimes(temp_var,eye(nx+nu));% Computing Sk

Q_k= blkdiag(R,Q); % Qk
h_k= -2*[R*ud;Q*xd]; %hk
c_k =(xd'*Q*xd) + (ud'*R*ud); % ck

if(k ==1)
    H_sum  = S_k'*Q_k*S_k;
    h_sum = h_k.'*S_k; 
    c_sum = c_k; 
else
    H_sum  = H_sum + S_k'*Q_k*S_k; 
    h_sum = h_sum +(h_k.'*S_k);
    c_sum = c_sum +c_k;
end
    
end

%%%% Handling the cost due to terminal penality %%%%
temp_Vec(1,N) =1;
S_N = otimes(temp_Vec,eye(nx+nu));

Q_N = blkdiag(R,Q); %Q-matrix due to terminal penality
h_N_temp= -2*[R*ud;(Q)*xd]; %h-matrix due to terminal penality
c_N = (xd'*(Q)*xd) + (ud'*R*ud); %c-matrix due to terminal penality

H_N = S_N'*Q_N*S_N;
h_N = h_N_temp'*S_N;


H = H_sum +H_N ;% Final H cost matrix
h = h_sum+h_N; % Final h cost matrix
c = c_sum + c_N; % Final c cost matrix
h=h'; %transpose the h to get the column vector

end