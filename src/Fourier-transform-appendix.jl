#_______________________________________________________________________________
#
# Code for Appendix: The Fourier Transform
#_______________________________________________________________________________

# VERSIONS: We used Julia Version 1.10.11, with the following package versions:
# [37e2e46d] LinearAlgebra
# [0c5d862f] Symbolics v6.46.0

# OVERVIEW: This document contains computations for verifying details in the
# appendix "The Fourier Transform", as well as code for converting between the
# probability coordinates and the Fourier q-coordinates that we use in the
# document. This document uses the same notation as the appendix. It is
# intended to be run line-by-line in a Julia repl session.

using LinearAlgebra, Symbolics
@variables p₀ p₃ p₂ p₁ p₄ # as defined in the appendix, these sum to 1.
p = [p₀, p₁, p₂, p₃, p₄]

# The specialized Fourier transform is the following matrix F. (Up to
# permutation of rows and columns, this is equivalent to what is provided in
# https://www.coloradocollege.edu/aapps/ldg/small-trees/small-trees_10.html).
F = [[1     1     1     1     1]
     [1     1  -1//3  -1//3  -1//3]
     [1  -1//3     1  -1//3  -1//3]
     [1  -1//3  -1//3     1  -1//3]
     [1  -1//3  -1//3  -1//3   1//3]]

# The following gives formulas for the Fourier q-coordinates in terms of p₀,
# p₁, p₂, p₃, p₄.
qAAA, qACC, qCAC, qCCA, qCGT = F*p
qAAA
qACC
qCAC
qCCA
qCGT

# We can invert the Fourier transformation as well, to obtain formulas for
# p₀,p₁, p₂, p₃, p₄ in terms of the q-coordinates:
@variables qAAA, qACC, qCAC, qCCA, qCGT
inv(F)*[qAAA, qACC, qCAC, qCCA, qCGT] # returns the vector p

# Define the transformation 𝓕
π_projection = [[0 1 0 0 0]
                [0 0 1 0 0]
                [0 0 0 1 0]
                [0 0 0 0 1]]
𝓕 = π_projection*F
𝓕*p # returns the simplified Fourier coordinates qACC qCAC, qCCA, qCGT

# Verification that the 4-dimensional Jacobian of T is sqrt(5):
@variables u₁ u₂ u₃ u₄
T = [1-u₁-u₂-u₃-u₄, u₁, u₂, u₃, u₄]
D = Symbolics.jacobian(T,[u₁,u₂,u₃,u₄])
det(D'*D) # The 4-dimensional Jacobian is the square root of this.


# Verification that the determinant of A is - 2*(4/3)^4:
A = (-4//3)*[[0 1 1 1]
            [1 0 1 1]
            [1 1 0 1]
            [1 1 1 1//2]]
det(A) - (-2*(4//3)^4) # Equals zero

# Finally, we verify that the restriction of 𝓕 to the simplex Δ₄ is the affine
# map 𝓕(x) = 1+Aπ(x) (equation 32 in the appendix). To verify this, observe
# that the following gives zero, since the probability coordinates sum to 1:
𝓕*p - ([1,1,1,1]+ A*π_projection*p)
