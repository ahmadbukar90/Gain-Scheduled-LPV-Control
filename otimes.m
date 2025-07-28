function M =otimes(A,B)

rA = size(A,1); % rows of matrix A
cA = size(A,2); % columns of matrix A

rB = size(B,1); % rows of matrix B
cB = size(B,2); % columns of matrix B
M = zeros (rA*rB,cA*cB); % Kronecker Matrix Dimensions

%% Initalizing the variables
var1 = 1;
var2 = 1;

for i=1 : size(A,1)
    for j = 1:size(A,2)
        temp_mat = A(i,j).*B;
        M(var1:i*rB,var2:j*cB) = temp_mat;
        var2 = var2+cB;
    end
    var1 = var1+rB;
    var2=1; %% rsetting the value of b
end



end
