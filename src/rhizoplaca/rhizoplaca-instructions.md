# Preliminaries

This file contains instructions for replicating the rhizoplaca example in the
paper. The example is based on the rhizoplaca network in Fig 4(a) in 

> Keuler, et al "Genome-scale data reveal the role of hybridization in
> lichen-forming fungi", Scientific Reports, (2020) 10:1497 |
> https://doi.org/10.1038/s41598-020-58279-x

The following instructions assume a Debian 12 operating system, but the steps
should be similar for MacOS and other Ubuntu/Debian-based Linux systems.

# Data preparation

First, obtain the gene tree dataset
[phylonet_BUSCO407_r4_haydenii.nex](../../data/phylonet_BUSCO407_r4_haydenii.nex)
by downloading the file `phylonet_BUSCO407_r4_haydenii.nex`
from
[https://doi.org/10.6084/m9.figshare.11299040](https://doi.org/10.6084/m9.figshare.11299040)
to `gitrepo/data`.

The source for this dataset is 

> Rachel Keuler, Alexis Garretson, Theresa Saunders, Robert Erickson, Nathan
> St.~Andre, Felix Grewe, H Thorsten Lumbsch, Jen-Pan Huang, Larry
> L.~St.~Claire, and Steven Leavitt. (2019). Genome-scale data reveal the role
> of hybridization in lichen-forming fungi Item. figshare. Dataset.
> https://doi.org/10.6084/m9.figshare.11299040.v2


We then edited the PhyloNet instruction block (located at the end of the nexus
file `phylonet_BUSCO407_r4_haydenii.nex`) in two ways:

- We changed the `maximum reticulations` parameter from 4 to 1 by replacing

  `InferNetwork_MPL (tree_1-tree_407) 4`
  with 
  `InferNetwork_MPL (tree_1-tree_407) 1`

- We excluded bootstrap sampling by removing the option flag `-b 50`.
  (Bootstrapping isn't necessary for our example, and including bootstrapping
  produced errors.)


# Instructions for replicating the network
To replicate the example network, do the following steps:

First, download PhyloNet version 3.8.2 (available at
[https://phylogenomics.rice.edu/html/phylonet.html](https://phylogenomics.rice.edu/html/phylonet.html)).
Save the file `PhyloNetv3_8_2.jar` to `gitrepo/src`.

Second, make sure java is installed correctly. For this example, I used the
folowing java version: `OpenJDK Runtime Environment (build
17.0.11+9-Debian-1deb12u1)`

Finally, run PhyloNet to estimate networks by running the following command
from the `/gitrepo`

```
java -jar src/PhyloNetv3_8_2.jar data/phylonet_BUSCO407_r4_haydenii.nex
```

This will call PhyloNet with maximum pseudolikelihood (as specified in the
nexus file), and return the top-five best scoring networks.

# Expected output

The best five networks we obtained by running the above code are

```
Inferred Network #1:
((((pari:1.0,(shus:1.0)#H1:1.0::0.6270893707082127):0.052149701388190514,((haydenii_complex:1.0,715f:1.0):0.42224547935742746,port:1.0):0.1601495353552912):0.4295167798034873,(mela:1.0,#H1:1.0::0.3729106292917873):0.1918552453196347):5.910309119954094,novo:1.0);
Total log probability: -3995168.2409823295
Visualize in Dendroscope : ((((pari,(shus)#H1),((haydenii_complex,715f),port)),(mela,#H1)),novo);

Inferred Network #2:
((((shus:1.0,(mela:1.0)#H1:1.0::0.5184465708981382):0.14832783102407038,(((715f:1.0,haydenii_complex:1.0):0.42325150252958577,port:1.0):0.16214828799388792,pari:1.0):0.1426851324490989):0.7001218931743756,#H1:1.0::0.48155342910186183):5.9115791163106675,novo:1.0);
Total log probability: -3995412.254519809
Visualize in Dendroscope : ((((shus,(mela)#H1),(((715f,haydenii_complex),port),pari)),#H1),novo);

Inferred Network #3:
((((shus:1.0)#H1:1.0::0.36399865381207996,mela:1.0):0.197218814627143,(((port:1.0,(haydenii_complex:1.0,715f:1.0):0.42353209308325496):0.159654043597711,#H1:1.0::0.63600134618792):0.0011774181844964955,pari:1.0):0.42956809268533375):5.910211696420865,novo:1.0);
Total log probability: -3995429.7294169627
Visualize in Dendroscope : ((((shus)#H1,mela),(((port,(haydenii_complex,715f)),#H1),pari)),novo);

Inferred Network #4:
(((mela:1.0,(shus:1.0)#H1:1.0::0.06847934236314923):5.933433279482921,(#H1:1.0::0.9315206576368508,(((715f:1.0,haydenii_complex:1.0):0.4234591024929822,port:1.0):0.16095369089584152,pari:1.0):0.12709066907043626):0.29358322025546174):5.9102129133179355,novo:1.0);
Total log probability: -3995565.6893627
Visualize in Dendroscope : (((mela,(shus)#H1),(#H1,(((715f,haydenii_complex),port),pari))),novo);

Inferred Network #5:
(((((((haydenii_complex:1.0,715f:1.0):0.3789656399904147,(port:1.0)#H1:1.0::0.9279940388410709):0.19785230712436316,pari:1.0):0.16718427000252375,shus:1.0):0.2753113253303149,mela:1.0):0.15025999127569925,#H1:1.0::0.0720059611589291):5.911380645369337,novo:1.0);
Total log probability: -3996584.3210581783
Visualize in Dendroscope : (((((((haydenii_complex,715f),(port)#H1),pari),shus),mela),#H1),novo);

```

There does not appear to be a seed parameter for this analysis, so the output
networks will vary slightly between runs. Nonetheless, the top-scoring network
we obtained (the first in the list) appears to closely match the network shown
in Figure 4(a) of [Keuler et al, 2020].

# Plotting the network

Instructions for plotting the top-two networks (shown in our paper) are found
in [rhizoplaca-plot.jl](rhizoplaca-plot.jl)

