--______________________________________________________________________________
--
-- This file contains Macaulay2 code for the computation in the proof of the
-- lemma in the appendix which shows the existence of a one-parameter real
-- family of solutions. It computes formulas for the 1-parameter fibers of the
-- parameterization map, in terms of the q-coordinates and the parameter a2.
--______________________________________________________________________________

-- Instructions: run this code with with Macaulay2 (version 1.25.06)

A = QQ[qACC, qCAC, qCCA, qCGT, a2]
K = frac(A);
B = K[a3, a4, a5, a6];
I = ideal(
      (              a2*a3*a5   ) - qACC,
      (   a3*a4*a5 +    a3*   a6) - qCAC,
      (a2*   a4    + a2*   a5*a6) - qCCA,
      (a2*a3*a4*a5 + a2*a3*a5*a6) - qCGT
    )
factor(a3 % I) -- formula for a3
factor(a4 % I) -- formula for a4
factor(a5 % I) -- formula for a5
factor(a6 % I) -- formula for a6
