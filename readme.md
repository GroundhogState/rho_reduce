# rho_reduce

Various algorithms for fast computation of partial trace operations, and benchmarks thereof.


## Getting Started

### Prerequisites

* MATLAB/[Octave](https://www.gnu.org/software/octave/) (Todo: Verify .m files work in octave)
* python 3.6

### Installing

Install with
```
git clone --recursive git://github.com/groundhogstate/rho_reduce.git
```
Update with
```
git submodule update --recursive --remote
```
## Running the tests

See [`TrC_examples.m`](lib/MATLAB/TrX_examples.m) linkfor examples in MATLAB syntax.

The alternative function `ptrace()` will allow the same syntax.

The two methods are compared in [`ptrace_vs_TrX.m`](https://github.com/GroundhogState/rho_reduce/blob/master/doc/fig/trace_time_matlab.png), with results

![](\doc\fig\trace_time_matlab.png)

## Usage



## About
See for Jonas Maziero's [paper](https://arxiv.org/abs/1609.00323) (\doc\mazerio_computing_partial_traces.pdf) which describes the implementation in ptrace.m.

The function TrX(rho,sys,dim) accepts the density matrix (rho) of a multipartite system with N elements each with dimension d_i, i=1:N. The total Hilbert space has therefore prod(d_i,i=1:N) dimensions. The function computes the trace 'over'  the systems specified by index in (sys), and returns a density matrix of the remaning systems ~(idx \cap sys). The vector (dim) specifies the dimensions, as required for the permutation of the density matrix into a product of non-square matrices, reducing the partial trace to partial inner product.

The function ptrace(rho,sys,dim) performs the same operation, but does so by computing the reduced state one matrix element at a time. Each element is produced in (I suspect) the optimal number of operations for promise-free well-conditioned matrices. It looks suspiciously vectorizable, and easily parallelized.

In the end, Toby's vectorized code beats my implementation of the iterative algorithm that Jonas describes. Toby's executes quickly in MATLAB because of the zero-cost commands reshape and permute, and calling LAPACK to compute the linear algebra. Todo: fix visual demo

However, there are many fewer operations in the iterative algorithm.

QETLAB has useful functions but may not be maintained any more. Could easily rewrite functions.
Currently, will just be useful for benchmarks in MATLAB.


## Built With

* PurpleBooth's [readme template](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* Nathaniel Johnston's [QETLAB](https://github.com/nathanieljohnston/QETLAB)
  * Tensor, ptrace
* Toby Cubitt's [TrX function](http://www.dr-qubit.org/matlab/TrX.m)

### Workflow
* [Atom](https://atom.io) extensions:
  * Hydrogen
  * markdown-preview-plus

## Contributing

Fork, edit, submit pull request. Don't break stuff.

## Versioning

## Authors

**Jacob Ross** [groundhogstate](https://github.com/groundhogstate)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

## TODO
* Port to python
* Port to C++
* Scaling analysis
* Test compiled versions
* Recursive many-body traces
* Benefits of parallelism?
