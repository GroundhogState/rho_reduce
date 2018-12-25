%The function TrX(rho,sys,dim) accepts the density matrix (rho) of a multipartite
% system with N elements each with dimension d_i, i=1:N. The total Hilbert
% space has therefore prod(d_i,i=1:N) dimensions. The function computes 
% the trace 'over'  the systems specified by index in (sys), and returns
% a density matrix of the remaning systems ~(idx \cap sys). The vector (dim)
% specifies the dimensions, as required for the permutation of the density
% matrix into a product of non-square matrices, reducing the partial trace
% to partial inner product. This executes quickly in MATLAB because of the
% zero-cost commands reshape and permute, and calling LAPACK to compute the
% linear algebra.

fprintf('=====Some examples of tracing=====\n')

% Tracing B out of qubits ABC returns a density matrix for the reduced
% states of AC.
% For a product of pure states:
A = rand_qubit();
B = rand_qubit();
C = rand_qubit();
ABC = Tensor(A,B,C);
fprintf('ABC is %u x %u, with trace %u\n',size(ABC),trace(ABC))
AC = TrX(ABC,2,[2,2,2])
% Trace out A and get density matrix for C back
C_red = TrX(AC,1,[2,2])

fprintf(['\n\nRetrieving the density matrix of a qubit from the product of\n',...
        'a 2 and 3-level system\n'])
D = rand_qubit();
E = toDM(normalize(exp(2j*pi*rand(3,1))));
DE = Tensor(D,E);
fprintf('DE is %u x %u with trace %u \n',size(DE),trace(DE))
TrX_E_DE= TrX(DE,2,[2,3])

% Retrieving the state of each of two three-level systems with a globally
% pure state
R = toDM(normalize(exp(2j*pi*randn(9,1))));
fprintf('R is %u x %u, with trace %u \n',size(R),trace(R))
P = TrX(R,2,[3,3]);
fprintf('P = Tr_Q (R) is %u x %u, with trace %u \n',size(P),trace(P))

