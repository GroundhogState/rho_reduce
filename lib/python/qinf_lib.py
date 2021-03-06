import numpy as np
import numpy.linalg as la
import warnings
import string

indices = string.ascii_lowercase

def dagger(op):
  return (op.transpose()).conj()

def toDM(psi):
    if len(np.shape(psi)) == 1: # accept only bra or ket
        if la.norm(psi) != 0:
            psi = psi/la.norm(psi) # normalize
        [row, col] = [psi[np.newaxis,:],psi[:,np.newaxis]];
        return np.dot(col,np.conj(row))


def Pauli(index):
    # Returns pauli I, X, Y, or Z corresponding to index (0 to 3)
    if index == 0:
        return np.eye(2) #Identity
    if index == 1:
        return np.array([[0,1],[1,0]]) # X
    if index == 2: # Y
        return np.array([[0,complex(0,-1)],[1,complex(0,1)]]) # Y
    if index == 3: # Z
        return np.array([[1,0],[0,-1]]) # Z
    else:
        warnings.warn('Invalid index to Pauli.py')


def qubit(P):
    # Returns density matrix of a qubit specified by polarization vector P
    if la.norm(P) != 0:
        P = P/la.norm(P)
    return 0.5*(Pauli(0) + sum([P[ii-1]*Pauli(ii) for ii in [1,2,3]]))

def partial_trace(state, n, modes):
    """
    Imported from https://github.com/XanaduAI/strawberryfields/blob/master/strawberryfields/backends/fockbackend/ops.py
    Computes the partial trace of a state over the modes in `modes`.
    Expects state to be in mixed state form.
    """
    left_str = [indices[2*i] + indices[2*i] if i in modes else indices[2*i:2*i+2] for i in range(n)]
    out_str = ['' if i in modes else indices[2*i:2*i+2] for i in range(n)] # Whence the twos?
    einstr = ''.join(left_str + ['->'] + out_str)

    return np.einsum(einstr, state)


def pytrace(rho,sys,dim):
    # A partial trace function.
    # Accepts one density matrix or state vector,  an array of dimensions of the subsystems,
    #       and an array of subsystems to trace out.
    #
    # Outputs the reduced density matrix obtained by tracing out the specified
    # subsystems.
    # Assumes:
    #   Input density matrix is Hermitian. Does this matter??
    # Depencies:
    #   None, although a tensor product function (like one from QETLab) is
    #   useful in applications.
    # References:
    #   Similar to Jonas Maziero's Fortran subroutine in https://github.com/jonasmaziero/LibForQ
    #       See also https://arxiv.org/abs/1601.07458
    #       Published in International Journal of Modern Physics CVol. 28, No. 01, 1750005 (2017)

    # Possible optimizations:
    #   - Trace out contiguous blocks internally if that buys any time
    #   - Convert to a general expression that reduces matrix in a single loop
    #   - Sort subsystems by size and trace out the large ones first
    #   - Replace MATLAB functions with smarter loops if more efficient
    #   - If you have absolute confidence in your inputs, remove the
    #       preconditioning steps to reduce total complexity.
    tr_rho = np.trace(rho)
    if tr_rho != 1:
        rho = rho/tr_rho                       #Normalizes density matrix
    N = len(dim)

    ## Function parameters
    tr_keep = set(np.r_[0:N]) - set(sys)       #The systems that will be KEPT - input is ones to trace out
    dim_out = np.prod(dim)                        #Will be reduced along the way
    left_keep = min(tr_keep)                   # The leftmost system that will remain. Used to condition the first stage.
    right_keep = max(tr_keep)                  # The rightmost system that will remain. Used to condition the first stage.
    dim_left = int(np.prod(dim[0:left_keep]))        # The dimension that will be traced out on the left
    dim_right = int(np.prod(dim[right_keep+1:]))    # The dimension that will be traced out on the right

    ## Function body
    # A small optimization: If there is a large contiguous block on either ,
    # trace over the largest first to reduce loop lens later.
    # Relegated to functions for readability & they're only called once.
    if (dim_left >= dim_right) & (left_keep > 0):
        # Trace left first
        dim_out = dim_out//dim_left
        rho_left = np.zeros([dim_out,dim_out])
        for i in np.r_[0:dim_out]:
            for j in np.r_[0:dim_out]:
                for p in np.r_[0:dim_left]:
                    rho_left[i,j] = rho_left[i,j]+rho[i+p*dim_out,j+p*dim_out]
        rho = rho_left
        if (right_keep<len(dim)):
            dim_out = dim_out//dim_right
            rho_right = np.zeros([dim_out,dim_out])
            for i in np.r_[0:dim_out]:
                for j in np.r_[0:dim_out]:
                    for p in np.r_[0:dim_right]:
                        rho_right[i,j]=rho_right[i,j]+rho[(i-1)*dim_right+p+1,(j-1)*dim_right+p+1]
            rho = rho_right
    elif (dim_right > dim_left) & (right_keep < len(dim)):
        # Trace right first
        dim_out = dim_out//dim_right
        rho_right = np.zeros([dim_out,dim_out])
        for i in np.r_[0:dim_out]:
            for j in np.r_[0:dim_out]:
                for p in np.r_[0:dim_right]:
                    rho_right[i,j]=rho_right[i,j]+rho[(i-1)*dim_right+p,(j-1)*dim_right+p]
        rho = rho_right
        if left_keep > 1:
            dim_out = dim_out//dim_left
            rho_left = np.zeros([dim_out,dim_out])
            for i in np.r_[0:dim_out]:
                for j in np.r_[0:dim_out]:
                    for p in np.r_[0:dim_left]:
                        rho_left[i,j] = rho_left[i,j]+rho[i+p*dim_out,j+p*dim_out]
            rho = rho_left
    # Performs traces over systems that aren't in contiguous blocks reaching to
    # the  of the system. Left as loop to reduce function call overheads.
    for k in np.r_[left_keep+1:right_keep]:
        if len(set([k])-set(sys)) == 0:
            d_mid = dim[k]
            dim_out = dim_out//d_mid
            d_low = np.prod(dim[k+1:right_keep+1])
            rho_mid = np.zeros([dim_out,dim_out])
            rho_mid = rho_mid.astype(complex)
            for i in np.r_[0:dim_out]:
                for j in np.r_[0:dim_out]:
                    ii = int(np.mod(i-1,d_low) + np.floor((i-1)/d_low)*d_mid*d_low)
                    jj = int(np.mod(j-1,d_low) + np.floor((j-1)/d_low)*d_mid*d_low)
                    for p in np.r_[0:d_mid]:
                        rho_mid[i,j] = rho_mid[i,j]+rho[1+ii+p*d_low,1+jj+p*d_low]
        rho = rho_mid
    return rho
