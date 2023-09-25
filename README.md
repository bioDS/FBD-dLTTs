# FBDsims
Simulations for identifiability paper

Castor was built using Rcpp - most underlying functionality is written in C and then imported into R. 
The simulations use a [forked version of castor](https://github.com/bioDS/castor), so that sampled ancestor trees can be produced.

To run the code using this version of castor, clone the repo linked above. 
Then, run `source updateRunFBD.sh` to build the package, and run the Rscript. You will need to make sure you do not have the standard castor implementation loaded into R.

If you have already built the package, you can run just the Rscript using `Rscript FBDensemble.R`.
