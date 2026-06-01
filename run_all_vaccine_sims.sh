#!/usr/bin/sh

scenario=A1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=A2 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=B1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=C1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=D1 start=1 Rscript ./scripts/vaccine_chunked.R
scenario=extra start=1 Rscript ./scripts/vaccine_chunked.R

