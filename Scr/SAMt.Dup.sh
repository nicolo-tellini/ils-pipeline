#!/bin/bash

BaseDir=$1

OutDir="MapDdp"
OutPath=$BaseDir"/"$OutDir
Nthreads=2

cd $BaseDir

if [[ ! -d $OutPath ]]; then
  mkdir -p $OutPath
fi

for Ind1 in $(ls MapRaw/*bam)
do
  Pref=$(echo $Ind1 | cut -d "/" -f 2 | cut -d "." -f 1,2)
  (
  samtools rmdup $Ind1 $OutPath/$Pref.srt.rmd.bam
  samtools index $OutPath/$Pref.srt.rmd.bam
  ) &
done

wait

