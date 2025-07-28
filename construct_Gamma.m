function [Phi,Psi] =construct_Gamma(A,B,N)

rA = size(A,1); %nx
cA = size(A,2); %nx

rB = size(B,1); %nx
cB = size(B,2); %nu

Gamma = zeros(N*rA,rA+N*(rA+cB)); %Constructing Gamma
In = -1*eye(rA,rA);

%% Developing the intermediate matrix
temp_mat = [A B In];
temp_r = size(temp_mat,1);
temp_c = size(temp_mat,2);

%% Initalizing the variables
a = 1;
b = 1;


for i = 1:N
    Gamma(a:i*temp_r,b:(b+temp_c-1)) = temp_mat;
    a = a+temp_r;
    b = b+(cA+cB);
end

Phi = Gamma(:,1:rA);
Psi = Gamma(:,rA+1:end);

end
