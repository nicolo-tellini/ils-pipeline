#!/bin/bash

Ref1Label=$1
BaseDir=$2

Nthreads=1
BamDir=${BaseDir}"/MapDdp"
DataExt=".rmd.bam"
OutDir=${BaseDir}"/VCFs/"
Nprocesses=20

cd $BamDir

if [[ ! -d $OutDir ]]; then
  mkdir -p $OutDir
else
  rm -rf $OutDir
  mkdir -p $OutDir
fi


for Ind1 in $(ls *${DataExt}) 
do
if (( i % Nprocesses == 0 )); then
 wait
fi
  ((i++))
 (

  (
  RefPath="/home/ntellini/pipeline-ILSWC/Ref/WE.genome.fa"
  RefName="WE"
  SampleName=$(echo $Ind1 | cut -d "." -f 1)

  samtools mpileup -u -min-MQ3 --output-tags AD,ADF,ADR,DP,SP --redo-BAQ -f $RefPath $Ind1 | bcftools call -vm -Oz > $OutDir"/"$SampleName"."$RefName".vcf.gz"
  wait
  tabix -p vcf $OutDir"/"$SampleName"."$RefName".vcf.gz" 
  )&
 )&
done
wait
