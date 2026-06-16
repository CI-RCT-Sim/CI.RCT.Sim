#!/usr/bin/sh

for i in {1..15}; do
  scenario=A1 row=$i Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..10}; do
  scenario=A2 row=$i Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..20}; do
  scenario=B1 row=$i Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..30}; do
  scenario=C1 row=$i Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..40}; do
  scenario=D1 row=$i Rscript ./scripts/vaccine_rowwise.R
done

for i in {1..45}; do
  scenario=extra row=$i Rscript ./scripts/vaccine_rowwise.R
done
