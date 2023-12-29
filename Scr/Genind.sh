#!/bin/bash

BaseDir=$1

AdmixtoolsVCFDir=${BaseDir}"/AdmixtoolsVCF/"

cd ${AdmixtoolsVCFDir}

for i in $(ls *.h.vcf.gz)
do
   indsamp=$(echo $i | cut -d"." -f1)

   cp -r ${BaseDir}"/Rep/WE_samp_SpEU_Sj.h.ind" .
   
   mv WE_samp_SpEU_Sj.h.ind ${indsamp}".h.ind"

   sed -i "s/\<Control\>/$indsamp/g" ${indsamp}".h.ind"

done

wait
