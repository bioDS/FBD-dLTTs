cd "/home/ket581/FBD_sims/FBD_version/castor"
Rscript ../updateAttributes.R
cd ../
R CMD build castor
R CMD INSTALL castor_1.7.11.tar.gz
Rscript FBDensemble.R
