# FBD-dLTTs
Determinstic lineage-through-time calculations for "The Fossilised Birth-Death Model is Identifiable"

The Figures folder contains the figures used in the paper, while the Graphing folder contains the R code and data used to reformat the output from running the FBDensemble.R file.

Castor was built using Rcpp - most underlying functionality is written in C and then imported into R. 
The simulations use a [forked version of castor](https://github.com/bioDS/castor), so that sampled ancestor trees can be produced.

Note that further work is necessary to make deterministic density calculations suitable for sampled ancestor trees.

To run the code using this version of castor, clone the repo linked above. 
Then, run `source updateRunFBD.sh` to build the package, and run the Rscript. You will need to make sure you do not have the standard castor implementation loaded into R.

If you have already built the package, you can run just the Rscript using `Rscript FBDensemble.R`.

