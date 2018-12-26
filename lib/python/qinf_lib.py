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
