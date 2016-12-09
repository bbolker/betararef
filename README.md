betararef
=========

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.162251.svg)](https://doi.org/10.5281/zenodo.162251)

Compute rarefied beta diversity and simulate communities (R package).  Accompanies Stier et al., "Isolating the effects of patch size and sampling effects on beta diversity via rarefaction", in prep.

**Installation**: for now, use `library(devtools); install_github("bbolker/betararef")` (don't forget `library("betararef")` to load the package). Then `help(package="betararef")` to see what's available.

`inst/batchfiles` contains example batch files used on a high-performance cluster (SHARCnet) to generate expected (realized) beta diversity across a broad range of community structures, sample sizes, and diversity indices.
