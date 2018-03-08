# Fast Atomic Norm Soft Thresholding (FastAST)
#### A fast primal-dual interior point method for line spectral estimation via atomic norm soft thresholding.

Implements the method of [1] for line spectral estimation via atomic norm minimization. If you use
this code, please cite this work.

[1] T. L. Hansen and T. L. Jensen, "A Fast Interior Point Method for Atomic Norm Soft Thresholding," submitted to *IEEE Transactions on Signal Processing*, 2018,
[preprint available on arXiv](https://arxiv.org/abs/1803.00862).

Abstract:
> The atomic norm provides a generalization of the l\_1-norm to continuous
> parameter spaces. When applied as a sparse regularizer for line spectral
> estimation the solution can be obtained by solving a convex optimization
> problem. This problem is known as atomic norm soft thresholding (AST). It can
> be cast as a semidefinite program and solved by standard methods. In the
> semidefinite formulation there are O(N^2) dual variables and a standard
> primal-dual interior point method requires at least O(N^6) flops per iteration.
> That has lead researchers to consider the alternating direction method of
> multipliers (ADMM) for the solution of AST, but this method is still somewhat
> slow for large problem sizes. To obtain a faster algorithm we reformulate AST
> as a non-symmetric conic program. That has two properties of key importance to
> its numerical solution: the conic formulation has only O(N) dual variables and
> the Toeplitz structure inherent to AST is preserved. Based on it we derive
> FastAST which is a primal-dual interior point method for solving AST. Two
> variants are considered with the fastest one requiring only O(N^2) flops per
> iteration. Extensive numerical experiments demonstrate that both variants of
> FastAST solve AST significantly faster than a state-of-the-art solver based on
> ADMM.


## Setup & Usage
A significant speedup of the code is obtained by building a mex version of the
generalized Schur algorithm which is used internally. Do so by running
`buildmex` in the MATLAB prompt. The MATLAB codegen feature is used to generate
mex files.

To allow static memory allocation, the largest N = length(y) that the generated
mex files will support is specified in `buildmex.m`.

MATLAB did not support recursion in codegen prior to version 9.0 (R2016a). On
these earlier versions a slower fallback approach is used which only uses
mex for the innermost iteration of the generalized Schur algorithm.

A simple example is provided in `example.m`.  Type `help solve_with_fastast` in
the MATLAB prompt for usage details.

This package redistributes the ADMM solver from
https://github.com/badrinarayan/astlinespec.

