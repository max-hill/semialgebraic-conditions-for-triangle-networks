#_______________________________________________________________________________
#
# This code provides the computation that the cutoff for distinguishability of
# rhizoplaca network is around s=.000177.
# _______________________________________________________________________________

# VERSIONS: This was run with Julia Version 1.10.11.

# INSTRUCTIONS: run the following code from a Julia repl

# functions to test model inclusion
function inM12p(a1,a2,a3,a4,a5,a6,δ)
    "returns true if the network is in the intersection of M1 and M2."
    N = a5*(δ*a4*(1-a3*a5) + a6*(1-δ)*(1-a3))
    D = (δ*a4*a5 + a6*(1-δ)) * (δ*a4*(1-a3*a5) + a5*a6*(1-δ)*(1-a3))
    return a1<N/D
end


function inM13p(a1,a2,a3,a4,a5,a6,δ)
    "returns true if the network is in the intersection of M1 and M3."
    N = a5*(δ*a4*(1-a2) + a6*(1-δ)*(1-a2*a5))
    D = (δ*a4 + a5*a6*(1-δ)) * (δ*a4*a5*(1-a2) + a6*(1-δ)*(1-a2*a5))
    return a1<N/D
end

# the estimated numerical parameters for the Rhizoplaca network:
δ = .373
s = .000177 # or .000178
ϵ = .0000001 # assume these to be of negligible length
d1, d2, d3, d4, d5, d6 = s, s, s, ϵ, .000192+.000482, ϵ
a1,a2,a3,a4,a5,a6 = exp.(-4/3*[d1,d2,d3,d4,d5,d6]) # convert to Fourier parameters

# check model inclusion
inM12p(a1,a2,a3,a4,a5,a6,δ)
inM13p(a1,a2,a3,a4,a5,a6,δ) 

# for s=.000177, the distribution is in M1 only. Rerunning the above with
# s=.000178, it is in the intersection of M1 and M3, and hence the hybrid node
# is not distinguishable.
