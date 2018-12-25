%The function ptrace(rho,sys,dim) accepts the density matrix (rho) of a multipartite
% system with N elements each with dimension d_i, i=1:N. The total Hilbert
% space has therefore prod(d_i,i=1:N) dimensions. The function computes 
% the trace 'over'  the systems specified by index in (sys), and returns
% a density matrix of the remaning systems ~(idx \cap sys). The vector (dim)
% specifies the dimensions, as required for the permutation of the density
% matrix into a product of non-square matrices, reducing the partial trace
% to partial inner product. This executes quickly in MATLAB because of the
% zero-cost commands reshape and permute, and calling LAPACK to compute the
% linear algebra.

fprintf('======== Some examples of tracing  with ptrace ========\n')



fprintf('\n==============================================\n')
fprintf('Retrieving the reduced states of each of two qubits in an entangled state')
Psi = [1 1/sqrt(2) 0 1/sqrt(2)]; % toDM normalizes automatically 
Rho = toDM(Psi/norm(Psi))
fprintf('\n\nThe reduced states are:\n')
psi1 = ptrace(Rho,2,[2,2])
psi2 = ptrace(Rho,1,[2,2])


fprintf('\n==============================================\n')
fprintf(['Tracing B out of a product of pure states ABC returns\n',...,
    '   a density matrix for the reduced states of AC\n'])
% For a product of pure states:
A = rand_qubit();
B = rand_qubit();
C = rand_qubit();
ABC = Tensor(A,B,C);

fprintf('ABC is %u x %u, with trace %u \n',size(ABC),trace(ABC))
AC = ptrace(ABC,2,[2,2,2])
fprintf('Trace out A and get density matrix for C back')
C_reduced = ptrace(AC,1,[2,2])
fprintf('\n ||C_reduced-C|| = %e \n',norm(C_reduced-C))





