#!/usr/bin/sh

echo installing renv
Rscript -e "install.packages('renv')"
echo restoring renv snapshot
Rscript -e "renv::restore()"
echo installing package from directory
Rscript -e "devtools::install()"
