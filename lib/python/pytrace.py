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
    tr_keep = set(np.r_[1:N+1]) - set(sys)       #The systems that will be KEPT - input is ones to trace out
    dim_out = np.prod(dim)                        #Will be reduced along the way
    left_keep = min(tr_keep)                   # The leftmost system that will remain. Used to condition the first stage.
    right_keep = max(tr_keep)                  # The rightmost system that will remain. Used to condition the first stage.
    dim_left = np.prod(dim[1:left_keep])        # The dimension that will be traced out on the left
    dim_right = np.prod(dim[right_keep+1:])    # The dimension that will be traced out on the right

    ## Function body
    # A small optimization: If there is a large contiguous block on either ,
    # trace over the largest first to reduce loop lens later.
    # Relegated to functions for readability & they're only called once.
    if (dim_left >= dim_right & left_keep > 1):
        # Trace left first
        dim_out = dim_out/dim_left
        rho_left = np.zeros(dim_out,dim_out)
        for i in np.r_(0,dim_out):
            for j in np.r_(0,dim_out):
                for p in np.r_(0,dim_left):
                    rho_left[i,j] = rho_left[i,j]+rho[i+p*dim_out,j+p*dim_out]
        rho = rho_left
        if (right_keep<len(dim)):
            #            print('hi')
            dim_out = dim_out/dim_right
            rho_right = np.zeros(dim_out,dim_out)
            for i in np.r_[0:dim_out]:
                for j in np.r_[0:dim_out]:
                    for p in np.r_[0:dim_right]:
                        rho_right[i,j]=rho_right[i,j]+rho[(i-1)*dim_right+p+1,(j-1)*dim_right+p+1]
            rho = rho_right
    elif dim_right > dim_left & right_keep < len(dim):
        # Trace right first
        dim_out = dim_out/dim_right
        rho_right = np.zeros(dim_out,dim_out)
        for i in np.r_[0:dim_out]:
            for j in np.r_[0:dim_out]:
                for p in np.r_[0:dim_right]:
                    rho_right[i,j]=rho_right[i,j]+rho[(i-1)*dim_right+p+1,(j-1)*dim_right+p+1]
        rho = rho_right
        if left_keep > 1:
            dim_out = dim_out/dim_left
            rho_left = np.zeros(dim_out,dim_out)
            for i in np.r_(0,dim_out):
                for j in np.r_(0,dim_out):
                    for p in np.r_(0,dim_left):
                        rho_left[i,j] = rho_left[i,j]+rho[i+p*dim_out,j+p*dim_out]
            rho = rho_left
    # Performs traces over systems that aren't in contiguous blocks reaching to
    # the  of the system. Left as loop to reduce function call overheads.
    for k in np.r_[left_keep+1:right_keep]:
        if ismember(k,sys):
            d_mid = dim[k]
            dim_out = dim_out/d_mid
            d_low = np.prod(dim[k+1:right_keep+1])
            rho_mid = zeros(dim_out)
            for i in np.r_[0:dim_out]:
                for j in np.r_[0:dim_out]:
                    ii = mod(i-1,d_low) + floor((i-1)/d_low)*d_mid*d_low
                    jj = mod(j-1,d_low) + floor((j-1)/d_low)*d_mid*d_low
                    for p in np.r_[0:d_mid]:
                        rho_mid[i,j] = rho_mid[i,j]+rho[1+ii+p*d_low,1+jj+p*d_low]
        rho = rho_mid
    return rho
