Fix of PENLAB bug (v.1.04)
==========================

Unfortunately, the latest version of PENLAB has a small bug due (probably) to MATLAB
changing the type of vector (row/column) returned by intersect/diff.

Replace the file `eval_alddx.m` in `source/@penlab` by the file provided in this folder.
