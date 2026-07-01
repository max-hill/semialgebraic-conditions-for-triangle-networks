#_______________________________________________________________________________
#
# This file contains code to make plots of the top-scoring Rhizoplaca network
#_______________________________________________________________________________

# VERSIONS: We used Julia Version 1.10.11, with the following package versions:
# [33ad39ac] PhyloNetworks v1.3.1
# [c0d5b6db] PhyloPlots v2.1.0
# [92933f4c] ProgressMeter v1.11.0
# [6f49c342] RCall v0.14.13
# R version 4.5.2 (platform: x86_64-pc-linux-gnu)

# INSTRUCTIONS: run this code from a Julia REPL in the directory
# /gitrepo/src/rhizoplaca/. The figures will be output to /gitrepo/figures/

using PhyloNetworks, PhyloPlots, RCall

# First, read in the highest-scoring network (the newick format was copied
# from instructions.md). In this case, the terminal branches and reticulation
# edges cannot be estimated from the datatype used, and by default PhyloNet
# assigns them value '1'. These lengths are not meaningful, so we delete them:
best_network = readnewick("((((pari,(shus)#H1:::0.6270893707082127):0.052149701388190514,((haydenii,715f):0.42224547935742746,port):0.1601495353552912):0.4295167798034873,(mela,#H1:::0.3729106292917873):0.1918552453196347):5.910309119954094,novo);")
rotate!(best_network,-4)
rotate!(best_network,-9)
plot(best_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)


# save the full network to a file
R"pdf('../../figures/top-rhizoplaca-network--coalescent-units.pdf', width=6, height=4)"
R"par(
  mar = c(2, 2, 2, 2),   # bottom, left, top, right (lines)
  oma = c(0, 0, 0, 0)    # outer margins
)"
plot(best_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)
R"mtext"("MPL inferred network (logL = -3995168)", line=-1)
R"dev.off()"


# next, restrict the network to only mela, shus, and pari:
three_leaf_network = readnewick("((((pari,(shus)#H1:::0.6270893707082127):0.052149701388190514,((haydenii_complex,715f):0.42224547935742746,port):0.1601495353552912):0.4295167798034873,(mela,#H1:::0.3729106292917873):0.1918552453196347):5.910309119954094,novo);")
deleteleaf!(three_leaf_network,8)
deleteleaf!(three_leaf_network,4)
deleteleaf!(three_leaf_network,5)
deleteleaf!(three_leaf_network,6)
rotate!(three_leaf_network,-9)
rotate!(three_leaf_network,-3)
plot(three_leaf_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)


# save the 3-leaf network to a file
R"pdf('../../figures/rhizoplaca-three-leaf--coalescent-units.pdf', width=6, height=4)"
R"par(
  mar = c(2, 2, 2, 2),   # bottom, left, top, right (lines)
  oma = c(0, 0, 0, 0)    # outer margins
)"
plot(three_leaf_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)
R"dev.off()"


# Next, we plot our networks with branch lengths in expected mutations per
# site. For conversion, we used the nucleotide diversity estimate
# π=2.5*10^(-3). This was obtained from Table 5 in Leavitt et al, J. Biogeogr.
# 2013, by averaging the estimates for π for the North American lineages of
# Rhizoplaca melanophthalma, C4b and C4d.
#
# This estimate allows us to convert branch lengths in coalescent units to
# expected mutations per site by multiplying by π/2≈10^-3 (the fungus are
# haploid). The result is the following network:
best_network = readnewick("((((pari,(shus)#H1:::0.0006270893707082127):0.000052149701388190514,((haydenii,715f):0.00042224547935742746,port):0.0001601495353552912):0.0004295167798034873,(mela,#H1:::0.0003729106292917873):0.0001918552453196347):5.000910309119954094,novo);")
rotate!(best_network,-4)
rotate!(best_network,-9)
plot(best_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)


# save the full network to a file
R"pdf('../../figures/top-rhizoplaca-network--evol-distance.pdf', width=6, height=4)"
R"par(
  mar = c(2, 2, 2, 2),   # bottom, left, top, right (lines)
  oma = c(0, 0, 0, 0)    # outer margins
)"
plot(best_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)
R"mtext"("MPL inferred network (logL = -3995168)", line=-1)
R"dev.off()"


# next, restrict the network to only mela, shus, and pari:
three_leaf_network = readnewick("((((pari,(shus)#H1:::0.0006270893707082127):0.000052149701388190514,((haydenii_complex,715f):0.00042224547935742746,port):0.0001601495353552912):0.0004295167798034873,(mela,#H1:::0.0003729106292917873):0.0001918552453196347):5.000910309119954094,novo);")
deleteleaf!(three_leaf_network,8)
deleteleaf!(three_leaf_network,4)
deleteleaf!(three_leaf_network,5)
deleteleaf!(three_leaf_network,6)
rotate!(three_leaf_network,-9)
rotate!(three_leaf_network,-3)
plot(three_leaf_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)


# save the 3-leaf network to a file
R"pdf('../../figures/rhizoplaca-three-leaf--evol-distance.pdf', width=6, height=4)"
R"par(
  mar = c(2, 2, 2, 2),   # bottom, left, top, right (lines)
  oma = c(0, 0, 0, 0)    # outer margins
)"
plot(three_leaf_network,showedgelength=true,showgamma=true,showedgenumber=false,shownodenumber=false)
R"dev.off()"
