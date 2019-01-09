# Correlation functions

A definitive application of the partial trace is the calculation of expectation values. For an operator $O$ and state $\rho$ we have

$$
\langle O\rangle = \textrm{Tr}(\rho O)
$$

Which is the projection of the state onto the operator subspace via the Hilbert-Schmidt inner product $\langle A,B\rangle = \textrm{Tr}(A^\dagger B)$, noting that $\rho$ and $O$ are both Hermitian. A full measurement is a projection onto a vector, and weaker measurements are projections onto higher-dimensional subspaces.

In particular, one may be interested in microscopic expectation values $\langle O_i\rangle = \textrm{Tr}(\rho O_i) =? \textrm{Tr}(\rho_i O_i)$,
intensive quantities $\langle\langle O\rangle\rangle = \langle \textrm{Tr}(O_i \rho)\rangle$, or n-point *connected* correlation functions $\langle O_i O_j\cdots O_n\rangle = \textrm{Tr} (\rho O_i O_j\cdots O_n)$, which one can compare with the  *disconnected* correlation functions $\prod \langle O_i\rangle $
