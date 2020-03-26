#!/bin/bash
ymd=`date +"%Y-%m-%d"`

animal_markers=( CO1 cox1 coxI COI )

for marker in "${animal_markers[@]}"
do
  echo ${marker}
  docker run -u $UID:$GID -v $PWD:/data --rm iimog/metabdb_dev\
   --outdir data/reference_library/uncurated/${ymd}/arthropoda_${marker}_full\
   --marker-search-string "${marker}"\
   --taxonomic-range Arthropoda\
   --sequence-length-filter 100:2000
 done
