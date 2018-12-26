# The partial trace

The *partial trace* is a central aspect of multilinear algebra with widespread physical applications. The intrinsic structure of the trace capture essential properties of quantum systems. The aim of this project is to understand the partial trace and its computation by implementing it in several ways, with the aim of finding the most effficient algorithmic implementation.

Todo: Probability perspective

## Preliminaries

Define two Gilbert spaces $F$ and $G$, with associated Banach spaces (of bounded Hermitian linear operators) $\mathcal{B}(F)$ and $\mathcal{B}(G)$. Define the tensor product of $F$ and $G$ by the map
$$
\otimes: F, G \to  F \otimes G =: H
$$

Which is linear in both arguments (other properties?). If $F$ and $G$ have dimension $d_F$ and $d_G$, then H has dimension $d_F\times d_G$. If $|f\rangle$ and $|g\rangle$ are orthnormal bases for $F$ and $G$ then the set of vectors $|f\rangle\otimes|g\rangle := |fg\rangle$ are an orthonormal basis for $H$.

The partial trace can then be written as a map
$$
Tr_F (O): F\otimes G \to G
$$

The density matrices of states in $F,G$ and $H$ are positive Hermitian operators with unit trace. The density matrices of $H$ can be seen as elements of $\mathbb{C}^{d_F d_G \times d_F d_G}$, leading to the definition of the partial trace of an operator $O$ as a partial inner product:


$$
Tr_F (O) = Tr_f O_{ij,kl} |f_i g_j\rangle\langle f_k g_l|
$$

$$
Tr_F (O) =  O_{ij,kl} \langle f_i | f_k \rangle\otimes |g_j\rangle\langle g_l|
$$

$$
Tr_F (O) =  O_{ij,kl} \delta_{ik}  |g_j\rangle\langle g_l|
$$

$$
Tr_F (O) =  O_{ij,il}  |g_j\rangle\langle g_l|
$$

Where $Tr_F (O)$ is an operator on $G$. The trace over $F$ is similar. The double-subscript on the matrix $H$ refers to the basis vector of $\mathcal{B}(H)$ formed by the product of the two corresponding basis vectors in $\mathcal{B}(F)$ and $\mathcal{B}(G)$.  

This notation suggests another alternative: Viewing operators as $(d_F \times d_F) \times (d_G \times d_G)$ tensors.


Which leads to the construction of the partial trace as an inner product of non-square matrices

$$
Tr_F(\rho) = \rho'^\dagger \rho'
$$
Where $\rho'$ is a $\frac{d_H}{d_G} \times d_G$ matrix, whose columns are indexed by the basis vectors of G and whose rows are indexed by the basis vectors spanning the space that will be traced out (this needs to be demystified in terms of conditional probabilities...)

Therefore to compute the partial trace of a bipartite system whose elements have dimension $d_F$ and $d_G$, use the algorithm:

```
rho_reshaped = reshape(rho,[d_H*d_H/d_G, d_G])
rho_G  = rho_reshaped'*rho_reshaped   % Note the conjugate transpose
```
Note: MATLAB reshapes by filling columns top-down, left-right, and Python reshapes by filling rows left-right, top-down. Therefore, the outputs would be transposes of one another, and the algorithm implemented accordingly. The above code is a MATLAB example.

Now, a subtlety: To use this method to retrieve the reduced state of $F$, it is a mistake to simply reshape as above with the dimensions interchanged: Indeed, for the case of two qubits, the result would be the same - $\rho_G!$.

The reason lies in the Hadamard-product structure of the basis for $H$. The reshaped density matrix must have the same indexing convention (col idx = basis vector to keep, col idx = vectors to trace out), and so the matrix rows and columns must be permuted to reflect this. In other words, one *can* use the prescription above for finding $\rho_F$, providing one swaps the positions of $F$ and $G!$

To see how, notice that the reshaping (defined above) breaks the operator into 'slices' of length $d_G$ - precisely the 'lowest-level' composition of the product basis, in the sense of the Hadamard product. Permuting the order of $F$ and $G$ is necessary so the reshaping procedure creates columns of `rho_reshaped` in 1-to-1 correspondence with distinguishable outcomes from measuring $F$.

One can simultaneously extend this to generality:

```
rho_factored = reshape(rho, (d_1,d_2,...d_N))
rho_permuted = permute(rho_factored(keep,traceout))
rho_reshaped = reshape(rho_permuted(dim_to_trace,dim_to_keep))
rho_reduced  = rho_reshaped'*rho_reshaped
```
## Conditional probabilities

QM is a theory of probabilities, though, so this discussion would be incomplete without the definition of the partial trace as an expectation value:
$$
Tr_F(O) = \sum_f p(g|f)p(f) O_g
$$




The first construction of the partial trace leads to the algorithm defined by J Maziero, and the second to the vectorized version implemented by T Cubitt. The Xanadu method, making use of `np.einsum()`, is a variation on the latter.
