#!/usr/bin/sh
rm -rf ~/.cache/R

echo installing renv
Rscript -e "install.packages('renv')"
Rscript -e "renv::restore(packages = 'renv')"
echo restoring renv snapshot
Rscript -e "renv::restore()"
echo installing devtools
Rscript -e "install.packages('devtools')"
echo installing package from directory
Rscript -e "devtools::install()"
echo run the sims
Rscript scripts/oncology_node4.R > first_try_msc.log 2>&1


