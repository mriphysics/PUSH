function [F] = optimisation_constraints(VOP,b_coeff,Nch,power_factor)

b_coeff = reshape(b_coeff,[],Nch);    
Nsp     = size(b_coeff,1);            
Nvop    = size(VOP,3);                

T = b_coeff' * b_coeff;

F = zeros(Nvop+Nch*Nsp,1);
num = 0;
for k = 1:Nvop
    num = num + 1;
    F(num) = real(power_factor * sum(sum(T .* VOP(:,:,k))));  
end

for k = 1:Nch*Nsp
    num = num + 1;
    F(num) = sqrt(real(b_coeff(k) * conj(b_coeff(k))));
end

end

