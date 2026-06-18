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

for i in {1..15}; do
  scenario=A1 row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..10}; do
  scenario=A2 row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..20}; do
  scenario=B1 row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..30}; do
  scenario=C1 row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..40}; do
  scenario=D1 row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..45}; do
  scenario=extra row=$i RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
done
