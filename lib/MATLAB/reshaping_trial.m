% rho_reshaped = reshape(rho,[d_H/d_G, d_G])
% rho_reduced  = rho_reduced'*rho_reduced   % Note the conjugate transpose
% psi = [1 1/sqrt(2) 1 1/sqrt(2)];
psi = [1 1/sqrt(2) 0 1/sqrt(2),1 0.5,0,-2j]; % toDM normalizes automatically 
rho = toDM(psi);
dim = [2 2 2];
D = prod(dim)
sys = 1;
N = length(dim);
sys_keep = setdiff(1:N,sys);
% rho_reshaped = reshape(rho,prod(dim)*prod(dim)/dim(2),dim(2))
% rho_reduced = rho_reshaped'*rho_reshaped
% TrX(rho,2,dim)
D = prod(dim);
rho_factored = reshape(rho,[dim,dim]);
rho_permuted = permute(rho_factored,[sys_keep,sys,sys_keep+N,sys+N]);
rho_reshaped = reshape(rho_permuted, (D^2)/prod(dim(sys_keep)),prod(dim(sys_keep)));
rho_reduced  = rho_reshaped'*rho_reshaped
TrX(rho,sys,dim)

% WORKS for two qubits. Breaks for three. Why? Something to do with the
% reshape/permute, the new pre-processing steps.