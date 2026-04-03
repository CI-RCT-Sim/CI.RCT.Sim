#!/usr/bin/sh

echo cloning git repo ===============================================
git clone -b sims_apr03 https://github.com/CI-RCT-Sim/CI.RCT.Sim.git CI.RCT.Sim
cd CI.RCT.Sim

echo installing renv ================================================
Rscript -e "install.packages('renv')"

echo restoring renv snapshot ========================================
Rscript -e "renv::restore()"

echo installing package from directory ==============================
Rscript -e "devtools::install()"
