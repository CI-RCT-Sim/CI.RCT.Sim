#!/usr/bin/sh

echo installing renv
Rscript -e "install.packages('renv')"
Rscript -e "renv::restore(packages = 'renv')"

echo restoring renv snapshot
Rscript -e "renv::restore()"

echo installing devtools
Rscript -e "install.packages('devtools')"

echo installing package from directory
Rscript -e "devtools::install()"

echo running row $row of scenario $scenario
RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R

