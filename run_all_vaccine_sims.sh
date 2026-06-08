#!/usr/bin/sh

# r script reads envrionment variables
# scenario A1, A2, ... scenario group
# start start at chunk 1, 2, ... to re-start sims from partial results

scenario=A1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=A2 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=B1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=C1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=D1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=extra start=1 Rscript ./scripts/vaccine_chunked.R

