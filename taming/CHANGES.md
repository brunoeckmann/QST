Changes
=======

September 27, 2017
------------------

- Started refactoring the former code for the (2,2,2) scenario into
  code that work on generic scenarios

- Added generation of random correlations
  
April 9, 2017
-------------

- Changed formulation in PENLAB Kullback-Leibler code to have probabilities
  always strictly positive (and avoid numerical issues when computing the KL
  divergence!)
  
- Added copy of Yanbao Zhang EM algorithm for tests

- Added PENLAB patch file to correct bug in v1.04

- Added test of consistency of KL/PENLAB with KL/EM.


April 8, 2017
-------------

- Added implementation of Kullback-Leibler minimization w.r.t. NNS or NPA sets
