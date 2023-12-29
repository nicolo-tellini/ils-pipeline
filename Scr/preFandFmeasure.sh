#!/bin/bash

BaseDir=$1

AdmixtoolsVCFDir=${BaseDir}"/AdmixtoolsVCF/"

cd ${AdmixtoolsVCFDir}

for i in $(ls *.h.vcf)
do
  indsampDir=$(echo $i | cut -d"." -f1)

  if [[ ! -d ${indsampDir} ]]; then
    mkdir -p ${indsampDir}
  fi

  cd ${indsampDir}

  wait

  #rm *

  wait

  cd ..

  cp -r ${i} ${indsampDir}
  
done

wait

# measure F
/usr/bin/time -v Rscript ${BaseDir}"/Scr/"f.R "${BaseDir}" >  ${BaseDir}"/Scr/Logs/f.out" 2>  ${BaseDir}"/Scr/Logs/Time.f.err" &

wait
