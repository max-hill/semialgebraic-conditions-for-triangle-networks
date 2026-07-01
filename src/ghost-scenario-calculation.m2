--_______________________________________________________________________________
--
-- This file contains Macaulay2 code for doing the computation in proof of
-- the ghost-lineage scenario theorem (It is doable by hand, but tedious).
--_______________________________________________________________________________

-- Macualay2 version 1.25.06

R = QQ[a1,a2,a3,a4,a5,a6,delta, MonomialOrder => Lex]
I = ideal(a6-a4*a5, a2-a1*a4)

-- In the proof of the theorem we showed that q∈M₁∩M₂ if and only if LHS<RHS,
-- where
LHS = a1*(delta*a4*a5+(1-delta)*a6)*(delta*a4*(1-a3*a5)+(1-delta)*a5*a6*(1-a3))
RHS = a5*(delta*a4*(1-a3*a5) + (1-delta)*a6*(1-a3))

-- We need to show that the following polynomial f is always negative:
f = (LHS-RHS)

-- Use the relations a6=a4*a5 and a2=a1*a4 to rewrite f without a1 and a4
f = factor (f % I)

-- Finally we show that f=h, where h is the (clearly negative) function
-- defined below:
h = -a6*(delta*(1-a2)*(1-a3*a5) + (1-delta)*a5*(1-a3)*(1-a2*a5))
f-h -- equals zero, so f=h. Since h is clearly negative, so is f.
