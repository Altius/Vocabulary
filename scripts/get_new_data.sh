#!/bin/bash

### Example usage:
# cat immunopedia-dnase-encode4-w-unactivated-0hrs.txt | cut -f 58 | tail -n +2 | ./get_new_data.sh - | gzip -c - > dat_immunopedia.txt.gz
# cat pancreatic-lines-production-qc-v2.txt | cut -f 58 | tail -n +2 | ./get_new_data.sh - | gzip -c - > dat_pancreatic.txt.gz

PEAKS_FILE=$1
INDEX="/home/meuleman/work/projects/ENCODE3/WM20180608_construct_masterlist_733samples/masterlist_DHSs_733samples_WM20180608_all.bed"

CMD="paste "
for PEAKS in `cat ${PEAKS_FILE}`; do
  CMD="$CMD <(bedmap --fraction-map 0.8 --indicator ${INDEX} ${PEAKS})"
done
eval $CMD


