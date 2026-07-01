#_______________________________________________________________________________
#
# This file contains code to replicate the plot in the section "Implications
# for Network Inference in Practice" and described in the final appendix.
# _______________________________________________________________________________

# VERSIONS: We used Julia Version 1.10.11, with the following package versions:
# [31c24e10] Distributions v0.25.125
# [b964fa9f] LaTeXStrings v1.4.0
# [91a5bcdd] Plots v1.41.6

# INSTRUCTIONS: run the code from a Julia repl in the directory 'gitrepo/src/'

using Plots, Plots.PlotMeasures, LaTeXStrings, Distributions
function inM1(q)
        (0<q[1]<1) &&
        (0<q[2]<1) &&
        (0<q[3]<1) &&
        (0<q[4]<1) &&
        (q[1]-q[4]>0) &&
        (q[4]-q[1]*q[2]>0) &&
        (q[1]*q[2]*q[3]-q[4]^2>0) &&
        (q[4]+q[1]*(q[4]-q[2]-q[3])>0)
end
inM12(q) = (q[4] + q[2]*(q[4]-q[1]-q[3])>0)
inM13(q) = (q[4] + q[3]*(q[4]-q[1]-q[2])>0)


function test_distinguishability_mc_a3(n,δ)
    "for a given δ, estimate the proportion of hybrid distinguishable networks
    for the N₁ rooted on edge a5, assuming a strict molecular clock"
    # initialize counters
    m12, m13, m123 = 0, 0, 0
    i=0
    while i<n
        # define random uniform time intervals h1, h2-h1, h3-h2, h4-h3
        h1 = rand(Uniform(0,.5))
        h2 = h1 + rand(Uniform(0,.5))
        h3 = h2 + rand(Uniform(0,.5))
        h4 = h3 + rand(Uniform(0,.5))
        # compute branch lengths (see figure 9 in the paper)
        t1 = h1
        t2 = h2
        t3 = h3
        t4 = h2 - h1
        t5 = 2*h4 - h2 - h3
        t6 = h3 - h1
        # compute Fourier parameters
        a1,a2,a3,a4,a5,a6 = exp.(-4/3*[t1,t2,t3,t4,t5,t6])
        # compute Fourier coordinates
        qACC = a2*a3*a5
        qCAC = δ*a1*a3   *a4*a5  +  (1-δ)*a1   *a3   *a6
        qCCA = δ*a1*a2   *a4     +  (1-δ)*a1*a2   *a5*a6
        qCGT = δ*a1*a2*a3*a4*a5  +  (1-δ)*a1*a2*a3*a5*a6
        q = [qACC, qCAC, qCCA, qCGT]
        # check if q lies in M1, M2, and M3.
        if inM1(q)
            i=i+1
            if inM12(q)
                m12=m12+1
            end
            if inM13(q)
                m13=m13+1
            end
            if inM12(q) && inM13(q)
                m123=m123+1
            end
        end
    end
    return 1-(m12 + m13 - m123)/n # return proportion distinguishable
end

# make the plot
x = (1:99)./100
y = [test_distinguishability_mc_a3(10000000, δ) for δ in x]
p = scatter(x,y,
            ylims=(0,.1), xlims=(0,1), xticks = 0:0.1:1, yticks = 0:0.01:1,
            size=(800,800), guidefontsize = 27, tickfontsize= 18,
            right_margin=3mm, bottom_margin=-3mm, legend=false,
            titlefontsize=27, markersize = 5, markershape = :circle,
            alpha = 1, xlabel = L"\mathrm{\delta\ (reticulation\ parameter)}",
            color = "blue", title = L"\mathrm{Proportion\ of\ samples\ in\ }\mathcal{M}_1\backslash(\mathcal{M}_2\cup\mathcal{M}_3)^c")

# save the plot
savefig(p,"../figures/distinguishability-a5-root-mc.pdf")
