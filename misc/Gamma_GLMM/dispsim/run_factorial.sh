#!/bin/bash
## Factorial over dispersion (shape = 1/phi) and cluster size, 20 clusters,
## RE sd 0.5, 100 replicates each, 12 workers. See README.md.
cd "$(dirname "$0")"
for nper in 4 20; do
    for shape in 5 1 0.5; do
        Rscript run_dispsim.R 100 20 $nper $shape 0.5 12 2>&1 | tail -1
        Rscript summarize_dispsim.R dispsim_g20_n${nper}_shape${shape}_sd0.5.rds \
            > dispsim_g20_n${nper}_shape${shape}_sd0.5_summary.txt
    done
done
