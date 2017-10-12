Matlab code to perform regularization of signaling distributions
================================================================

This code implements regularization of signaling distributions by a variety of methods.

Projection
----------

Provided by `FindProjection.m`, does not require any solver.

Nearest approximation with respect to KL divergence (PENLAB)
------------------------------------------------------------

Requires: [PENLAB](http://web.mat.bham.ac.uk/kocvara/penlab/) with the provided
patch (see `/fix_penlab_bug`).

The function `FindNA_KL_PENLAB.m` provides an implementation of
KL minimization with respect to:

- the nonsignaling polytope, represented by linear inequalities (level = []),
- semidefinite outer approximations of the quantum set (proper hierarchy levels),

using the nonlinear solver PENLAB.
