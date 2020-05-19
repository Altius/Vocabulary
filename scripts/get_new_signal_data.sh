#!/bin/bash

### Example usage:
# cat immunopedia-dnase-encode4-w-unactivated-0hrs.txt | cut -f 62 | tail -n +2 | ./get_new_signal_data.sh - | gzip -c - > dat_immunopedia_signal.txt.gz
# cat pancreatic-lines-production-qc-v2.txt | cut -f 63 | tail -n +2 | ./get_new_signal_data.sh - | gzip -c - > dat_pancreatic_signal.txt.gz

SIGNAL_FILE=$1
INDEX="/home/meuleman/work/projects/ENCODE3/WM20180608_construct_masterlist_733samples/masterlist_DHSs_733samples_WM20180608_all_indexIDs.txt"

CMD="paste "
for SIGNAL in `cat ${SIGNAL_FILE}`; do
  CMD="$CMD <(awk 'BEGIN{OFS=\"\t\"}{print \$1, \$9, \$9+1, \$4, \$5}' ${INDEX} | bedmap --max - ${SIGNAL})"
done
eval $CMD
#echo $CMD


