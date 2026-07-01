#_______________________________________________________________________________
#
# This file contains all code for the simulation in the Section "The Volume of
# Model Intersections" and code to replicate the plot of the figure.
# _______________________________________________________________________________

# VERSIONS: We used Julia Version 1.10.11, with the following package versions:
# [31c24e10] Distributions v0.25.125
# [b964fa9f] LaTeXStrings v1.4.0
# [91a5bcdd] Plots v1.41.6
# [37e2e46d] LinearAlgebra

# INSTRUCTIONS: run this code in a Julia REPL in the directory 'gitrep/src/'.
# Plots will be saved to 'gitrepo/figures/'


# Load packages and define utility functions
using Distributions, LinearAlgebra, Plots, Plots.PlotMeasures, LaTeXStrings

# function to test whether a point q is in M₁ 
function inM1(q)
    "return true if q ∈ M₁"
        (0<q[1]<1) &&
        (0<q[2]<1) &&
        (0<q[3]<1) &&
        (0<q[4]<1) &&
        (q[1]-q[4]>0) &&
        (q[4]-q[1]*q[2]>0) &&
        (q[1]*q[2]*q[3]-q[4]^2>0) &&
        (q[4]+q[1]*(q[4]-q[2]-q[3])>0)
end

function inM2(q)
    "return true if q ∈ M₂"
        (0<q[1]<1) &&
        (0<q[2]<1) &&
        (0<q[3]<1) &&
        (0<q[4]<1) &&
        (q[2]-q[4]>0) &&
        (q[4]-q[1]*q[2]>0) &&
        (q[1]*q[2]*q[3]-q[4]^2>0) &&
        (q[4] + q[2]*(q[4]-q[1]-q[3])>0)
end

function inM3(q)
    "return true if q ∈ M₃"
        (0<q[1]<1) &&
        (0<q[2]<1) &&
        (0<q[3]<1) &&
        (0<q[4]<1) &&
        (q[3]-q[4]>0) &&
        (q[4]-q[1]*q[3]>0) &&
        (q[1]*q[2]*q[3]-q[4]^2>0) &&
        (q[4]+q[3]*(q[4]-q[1]-q[2])>0)
end

# functions to test whether a point q ∈ M₁ is also in M₂ or M₃
inM12(q) = (q[4] + q[2]*(q[4]-q[1]-q[3])>0)
inM13(q) = (q[4] + q[3]*(q[4]-q[1]-q[2])>0)


#_______________________________________________________________________________
#
# Code for Observation 1 in the paper (model overlap in the simplex)
#_______________________________________________________________________________

# The following function was used to compute the values in the Venn diagram
# figure:

function estimateVolumes(n)
    "estimate the proportion of the simplex corresponding to each of the
    modems M₁,M₂, and M₃, and their intersections."
    # initialize counters
    m1,m2,m3,m12,m13,m23,m123=0,0,0,0,0,0,0
    # define the Fourier transform (this is 𝓕 in Appendix B)
    𝓕 = [[1     1  -1/3  -1/3  -1/3]
         [1  -1/3     1  -1/3  -1/3]
         [1  -1/3  -1/3     1  -1/3]
         [1  -1/3  -1/3  -1/3   1/3]]
    for i in 1:n
        # draw a point uniformly from 4d simplex
        p = rand(Dirichlet([1,1,1,1,1])) 
        # Take the Fourier transform, to obtain the vector q=[qACC,qCAC,qCCA,qCGT]
        q=𝓕*p
        # Count instances of model inclusion
        inM1(q)
        m1=m1+inM1(q)
        m2=m2+inM2(q)
        m3=m3+inM3(q)
        m12=m12+(inM1(q) && inM2(q))
        m13=m13+(inM1(q) && inM3(q))
        m23=m23+(inM2(q) && inM3(q))
        m123=m123+(inM1(q) && inM2(q) && inM3(q))
    end
    # estimate proportions
    pM1, pM2, pM3, pM12, pM13, pM23, pM123 = m1/n, m2/n, m3/n, m12/n, m13/n, m23/n, m123/n
    # print percentage results
    print("\n\nPROPORTIONS\n")
    println("Proportion in M₁          = ", round(pM1*100,digits=3),"%")
    println("Proportion in M₂          = ", round(pM2*100,digits=3),"%")
    println("Proportion in M₃          = ", round(pM3*100,digits=3),"%")
    print("\n")
    println("Proportion in M₁∩M₂       = ", round(pM12*100,digits=3),"%")
    println("Proportion in M₁∩M₃       = ", round(pM13*100,digits=3),"%")
    println("Proportion in M₂∩M₃       = ", round(pM23*100,digits=3),"%")
    print("\n")
    println("Proportion in M₁∩M₂∩M₃    = ", round(pM123*100,digits=3),"%")
    print("\n")
    println("Proportion in M₁∪M₂∪M₃    = ", round((pM1+pM2+pM3-pM12-pM13-pM23+pM123)*100,digits=3),"%")
    println("Proportion in (M₁∪M₂∪M₃)ᶜ = ", round((1-(pM1+pM2+pM3-pM12-pM13-pM23+pM123))*100,digits=3),"%")
    print("\n")
    println("Proportion in (M₁∩M₂)\\M₃  = ", round((pM12 - pM123)*100,digits=3),"%")
    println("Proportion in (M₁∩M₃)\\M₂  = ", round((pM13 - pM123)*100,digits=3),"%")
    println("Proportion in (M₂∩M₃)\\M₁  = ", round((pM23 -  pM123)*100,digits=3),"%")
    print("\n")
    println("Proportion in M₁\\(M₂∪M₃)  = ", round((pM1 - pM12 - pM13 + pM123)*100,digits=3),"%")
    println("Proportion in M₂\\(M₁∪M₃)  = ", round((pM2 - pM12 - pM23 + pM123)*100,digits=3),"%")
    println("Proportion in M₃\\(M₁∪M₂)  = ", round((pM3 - pM13 - pM23 + pM123)*100,digits=3),"%")
    print("\n")
    println("Proportion in M₁∩M₂\\M₃    = ", round((pM12 - pM123)*100,digits=3),"%")
    println("Proportion in M₁∩M₃\\M₂    = ", round((pM13 - pM123)*100,digits=3),"%")
    println("Proportion in M₂∩M₃\\M₁    = ", round((pM23 - pM123)*100,digits=3),"%")
    print("\n")
    println("Proportion of points in M₁ that are also in M₂∪M₃    = ",
            round((pM12 + pM13 - pM123)/pM1,digits=3),"%")
    # print volume results
    c = 2/sqrt(5)*(4/3)^4
    d = sqrt(5)/24
    print("\n\nVOLUMES\n")
    println("Volume of M₁          = ", round(pM1*c*d,digits=5))
    println("Volume of M₂          = ", round(pM2*c*d,digits=5))
    println("Volume of M₃          = ", round(pM3*c*d,digits=5))
    print("\n")
    println("Volume of M₁∩M₂       = ", round(pM12*c*d,digits=5))
    println("Volume of M₁∩M₃       = ", round(pM13*c*d,digits=5))
    println("Volume of M₂∩M₃       = ", round(pM23*c*d,digits=5))
    print("\n")
    println("Volume of M₁∩M₂∩M₃    = ", round(pM123*c*d,digits=5))
    print("\n")
    println("Volume of 𝓕(Δ₄)      = ", round(c*d,digits=5))
    println("Volume of M₁∪M₂∪M₃    = ", round((pM1+pM2+pM3-pM12-pM13-pM23+pM123)*c*d,digits=5))
    println("Volume of (M₁∪M₂∪M₃)ᶜ = ", round((1-(pM1+pM2+pM3-pM12-pM13-pM23+pM123))*c*d,digits=5))
    print("\n")
    println("Volume of (M₁∩M₂)\\M₃  = ", round((pM12 - pM123)*c*d,digits=5))
    println("Volume of (M₁∩M₃)\\M₂  = ", round((pM13 - pM123)*c*d,digits=5))
    println("Volume of (M₂∩M₃)\\M₁  = ", round((pM23 -  pM123)*c*d,digits=5))
    print("\n")
    println("Volume of M₁\\(M₂∪M₃)  = ", round((pM1 - pM12 - pM13 + pM123)*c*d,digits=5))
    println("Volume of M₂\\(M₁∪M₃)  = ", round((pM2 - pM12 - pM23 + pM123)*c*d,digits=5))
    println("Volume of M₃\\(M₁∪M₂)  = ", round((pM3 - pM13 - pM23 + pM123)*c*d,digits=5))
    print("\n")
    println("Volume of M₁∩M₂\\M₃    = ", round((pM12 - pM123)*c*d,digits=5))
    println("Volume of M₁∩M₃\\M₂    = ", round((pM13 - pM123)*c*d,digits=5))
    println("Volume of M₂∩M₃\\M₁    = ", round((pM23 - pM123)*c*d,digits=5))
end

estimateVolumes(1000000000) # may take a few minutes



#_______________________________________________________________________________
#
# Code for Observation 2 in the paper (model overlap in the parameter space)
#_______________________________________________________________________________


# The following function was used to generate networks of the form N₁, with
# parameters drawn uniformly at random:

function estimate_distinguishable_proportion_uniform_case(n)
    "Estimate the proportion of networks in M₁, M₁∩M₂, M₁∩M₃, and M₁∩M₂∩M₃
    when all numerical parameters are drawn unifromly at random from the
    parameter space."
    m12, m13, m123 = 0, 0, 0
    i=0
    while i<n
        # generate numerical parameters uniformly at random
        a1,a2,a3,a4,a5,a6,δ = rand(Uniform(0,1),7)
        # compute Fourier parameterization
        qACC = a2*a3*a5
        qCAC = δ*a1*a3   *a4*a5  +  (1-δ)*a1   *a3   *a6
        qCCA = δ*a1*a2   *a4     +  (1-δ)*a1*a2   *a5*a6
        qCGT = δ*a1*a2*a3*a4*a5  +  (1-δ)*a1*a2*a3*a5*a6
        q = [qACC, qCAC, qCCA, qCGT]
        # record model inclusions
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
    # print results
    println("Proportion in M₁∩M₂     = ", m12/n)
    println("Proportion in M₁∩M₃     = ", m13/n)
    println("Proportion in M₁∩M₂∩M₃  = ", m123/n)
    println("Prop. indistinguishable = ", (m12 + m13 - m123)/n)
    println("Prop. distinguishable   = ", 1-(m12 + m13 - m123)/n)    
end

estimate_distinguishable_proportion_uniform_case(1000000000)


# The following code was used to generate the multiple bell-curve plot

function estimate_distinguishable_proportion(n,δ,m)
    "Given a reticulation parameter δ, estimate the proportion of networks
    that satisfy both distinguishability criteria (i.e., are in M1 but not in
    M2 or M3). The estimate is based on n randomly generated networks with
    branch lengths drawn uniformly from (0,M)"
    m12, m13, m123 = 0, 0, 0
    i=0
    while i<n
        # generate branch lengths
        d1,d2,d3,d4,d5,d6 = rand(Uniform(0,m),6)
        a1,a2,a3,a4,a5,a6 = exp.(-4/3*[d1,d2,d3,d4,d5,d6])
        # compute Fourier parameterization
        qACC = a2*a3*a5
        qCAC = δ*a1*a3   *a4*a5  +  (1-δ)*a1   *a3   *a6
        qCCA = δ*a1*a2   *a4     +  (1-δ)*a1*a2   *a5*a6
        qCGT = δ*a1*a2*a3*a4*a5  +  (1-δ)*a1*a2*a3*a5*a6
        q = [qACC, qCAC, qCCA, qCGT]
        # record model inclusions
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
    # return proportion of networks in which the hybrid is distinguishable
    return 1-(m12 + m13 - m123)/n
end

# generate data for the plot
x = (1:99)./100
n= 10000000
y1 = [estimate_distinguishable_proportion(n, δ, 1) for δ in x]
y2 = [estimate_distinguishable_proportion(n, δ, .5) for δ in x]
y3 = [estimate_distinguishable_proportion(n, δ, .25) for δ in x]
y4 = [estimate_distinguishable_proportion(n, δ, .05) for δ in x]
y5 = [estimate_distinguishable_proportion(n, δ, .01) for δ in x]

# make the plot
p = scatter(x,y1,
            ylims=(0,.008), xlims=(0,1), xticks = 0:0.1:1, yticks = 0:0.001:1,
            size=(800,800), guidefontsize = 18, legendfontsize = 18,
            tickfontsize= 13, right_margin=3mm, bottom_margin=-3mm,
            titlefontsize=18, marker=:cross, markerstrokewidth=0.5,
            #xlabel = L"\mathrm{\delta\ (reticulation\ parameter)}",
            #ylabel = L"\mathrm{Proportion\ with\ distinguishable\ hybrid\ node}",
            label= L"m=1")
scatter!(p, x, y2, label = L"m=0.5", marker=:xcross, markerstrokewidth = 0.3)
scatter!(p, x, y3, label = L"m=0.25", marker=:star4, markerstrokewidth = 0.3)
scatter!(p, x, y4, label = L"m=0.05", marker=:diamond, markerstrokewidth = 0.3)
scatter!(p, x, y5, label = L"m=0.01", marker=:circle, markerstrokewidth = 0.3)

# save the plot
savefig(p,"../figures/plot-for-observation-2.pdf")

