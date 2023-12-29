#!/bin/bash

BaseDir=$1

/usr/bin/time -v Rscript "${BaseDir}"/Scr/GenChrDir.R "${BaseDir}" > "${BaseDir}"/Scr/Logs/GenChrDir.out 2> "${BaseDir}"/Scr/Logs/Time.GenChrDir.err &

wait

####### insert the header

AdmixtoolsVCFDir=${BaseDir}"/AdmixtoolsVCF/"

cd ${AdmixtoolsVCFDir}

for i in $(ls -d */ | grep -v "results")
do
   cd $i

  for j in $(ls -d */)
   do
  
   samp=$(echo $j | cut -d"_" -f1)

   chr=$(echo $j | cut -d"_" -f2 | cut -d"/" -f1)
   
   cd $j

   cp -r ${BaseDir}"/Rep/header.txt" .

   cat header.txt $samp"."$chr".vcf" > $samp"."$chr".h.vcf"

   rm $samp"."$chr".vcf"  header.txt

   cat $samp"."$chr".h.vcf" | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > out_sorted.vcf

   rm $samp"."$chr".h.vcf"

   mv out_sorted.vcf $samp"."$chr".h.vcf"

   cd ..    

   done

 cd ..
 
done

wait

# measure DF

/usr/bin/time -v Rscript "${BaseDir}"/Scr/measureDF.R "${BaseDir}" > "${BaseDir}"/Scr/Logs/measureDF.out 2> "${BaseDir}"/Scr/Logs/Time.measureDF.err &
