#!/bin/bash

BaseDir=$1
Nprocesses=20

SNPDir=${BaseDir}"/AdmixtoolsVCF/"

cd ${SNPDir}

for j in $(ls *.plink.vcf)
do
if (( i % Nprocesses == 0 )); then
 wait
fi
  ((i++))
(
samp=$(echo $j | cut -d"." -f1)

cat ${BaseDir}"/Rep/"header.txt ${SNPDir}"/"$j > ${SNPDir}"/"${samp}".h.vcf"
) &
done

wait
