#!/usr/bin/sh

# install packages, restore renv snapshot, etc.
echo installing renv
Rscript -e "install.packages('renv')"
echo restoring renv snapshot
Rscript -e "renv::restore()"
echo installing package from directory
Rscript -e "devtools::install()"

# run scenarios
# enviornment variables
#   * scenario: scenario class
#   * SLURM_ARRAY_TASK_ID: row-number in scenario class (automatically set when submitting as slurm job array, manually set here)
#   * RENV_CONFIG_SANDBOX_ENABLED: disable renv library sandbox (speeds up initialisation of parallel cluster by avoiding simultaneous file system reads)
#
# rows in scenario classes:
# A1 1-15
# A2 1-10
# B1 1-20
# C1 1-30
# D1 1-40
# extra 1-45
scenario=A1    SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
scenario=A2    SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
scenario=B1    SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
scenario=C1    SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
scenario=D1    SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
scenario=extra SLURM_ARRAY_TASK_ID=1 RENV_CONFIG_SANDBOX_ENABLED=FALSE Rscript ./scripts/vaccine_rowwise.R
