#!/usr/bin/env bash

set -euo pipefail

IN="../data/datasets/hca_brain-organoids_grouped_restricted/coverages_chrombpnet/ualf_Glia_unstranded.bw"
OUT="../data/processed-data/hca_brain-organoids_grouped_restricted/chromatin-accessibility/bins_500/chrombpnet/ualf_Glia.bw"
BWTOOLS="${HOME}/apps/make/bwtool/bwtool"

mkdir -p $(dirname "$OUT")
if [ -e $OUT ]; then
	rm $OUT
fi

$BWTOOLS summary 500 $IN $OUT -with-sum -fill=0
