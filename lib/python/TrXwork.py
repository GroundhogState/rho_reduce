import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import qinf_lib as qi
import string

indices = string.ascii_lowercase

def ptr(state,n,modes):
    # n = number of modes (subsystems in qubit case?)
    left_str = [indices[2*i] + indices[2*i] if i in modes else indices[2*i:2*i+2] for i in range(n)]
    out_str = ['' if i in modes else indices[2*i:2*i+2] for i in range(n)] # Whence the twos?
    einstr = ''.join(left_str + ['->'] + out_str)
    print(einstr)
    # return left_str
    print(np.shape(state))
    return np.einsum(einstr, state)

print('====== Some examples of pythonic tracing')
print('== Retrieving the reduced states of an entangled qubit pair:')
print('==   Reshape-and-sum method.')
psi = np.array([1,1/np.sqrt(2),0,1/np.sqrt(2)])
rho = qi.toDM(psi)
dim = [2,2]
# Reshape into the form dim_keep, dim_traceout
# Trace out the SECOND system and get the FIRST back
rho_reshaped = np.reshape(rho,(dim[1],(np.prod(dim)**2)//dim[1]))
rho_reduced = np.dot(rho_reshaped,qi.dagger(rho_reshaped))
print(rho_reduced)




"""
Test case
    psi = np.array([1,0,1,1])
    rho = qi.toDM(psi)
    ptr(rho,3,[1]) produces
    einstr = abccef-> abef
    ptr(rho,2,[1]) produces
    einstr = abcc->ab   <= Edge case bug?!
"""

### using np.einsum
# a = np.arange(25).reshape(5,5)
#b = np.arange(5)
#c = np.arange(6).reshape(2,3)
#d=np.einsum('ii', a)
#print(a)
#print(d)
