# Preliminaries

This file contains instructions for replicating the rhizoplaca example in the
paper. The example is based on the rhizoplaca network in Fig 4(a) in 

> Keuler, et al "Genome-scale data reveal the role of hybridization in
> lichen-forming fungi", Scientific Reports, (2020) 10:1497 |
> https://doi.org/10.1038/s41598-020-58279-x

The following instructions assume a Debian 12 operating system, but the steps
should be similar for MacOS and other Ubuntu/Debian-based Linux systems.

# Data preparation
To obtain and prepare the dataset, there are two steps:

First, obtain the gene tree dataset `phylonet_BUSCO407_r4_haydenii.nex` by
downloading it from
[https://doi.org/10.6084/m9.figshare.11299040](https://doi.org/10.6084/m9.figshare.11299040),
and saving it to `gitrepo/data`.

The source for this dataset is 

> Rachel Keuler, Alexis Garretson, Theresa Saunders, Robert Erickson, Nathan
> St. Andre, Felix Grewe, H Thorsten Lumbsch, Jen-Pan Huang, Larry
> L. St. Claire, and Steven Leavitt. (2019). Genome-scale data reveal the role
> of hybridization in lichen-forming fungi Item. figshare. Dataset.
> https://doi.org/10.6084/m9.figshare.11299040.v2


Second, edit the PhyloNet instruction block (located at the end of the nexus
file `phylonet_BUSCO407_r4_haydenii.nex`) in two ways:

- Change the `maximum reticulations` parameter from 4 to 1 by replacing

  `InferNetwork_MPL (tree_1-tree_407) 4`
  with 
  `InferNetwork_MPL (tree_1-tree_407) 1`

- Exclude bootstrap sampling by removing the option flag `-b 50`.
  (Bootstrapping isn't necessary for our example, and including
  bootstrapping produced errors.)


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

The top-scoring five networks we obtained by running the above code are in [rhizoplaca-best-networks.txt](../../data/rhizoplaca-best-networks.txt). Our best network was

```
Inferred Network #1:
((((pari:1.0,(shus:1.0)#H1:1.0::0.6270893707082127):0.052149701388190514,((haydenii_complex:1.0,715f:1.0):0.42224547935742746,port:1.0):0.1601495353552912):0.4295167798034873,(mela:1.0,#H1:1.0::0.3729106292917873):0.1918552453196347):5.910309119954094,novo:1.0);
Total log probability: -3995168.2409823295
Visualize in Dendroscope : ((((pari,(shus)#H1),((haydenii_complex,715f),port)),(mela,#H1)),novo);
```

There does not appear to be a seed parameter for this analysis, so the output
networks will vary slightly between runs. Nonetheless, the top-scoring network
we obtained (the first in the list) appears to closely match the network shown
in Figure 4(a) of [Keuler et al, 2020].

# Plotting the network

Instructions for plotting the top-two networks (shown in our paper) are found
in [rhizoplaca-plot.jl](rhizoplaca-plot.jl)

