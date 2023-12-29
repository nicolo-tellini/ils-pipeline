#!/bin/bash

#####################
### user settings ###
#####################
 
# assemblies labels (must match prefix of assemblies stored in "Rep" folder)
Ref1Label="WE"
Ref1="WE"

#####################
### settings' end ###
#####################

### part 0

# base folder
BaseDir=$(dirname "$(pwd)")

# check Logs folder
if [[ ! -d Logs ]]; then mkdir Logs; fi

### part I: references & annotations initialization

#mkdir -p "${BaseDir}"/Ref/Ann
#mkdir -p "${BaseDir}"/Ref/Asm
#cp "${BaseDir}"/Rep/Asm/"${Ref1Label}".genome.fa "${BaseDir}"/Ref &
#cp "${BaseDir}"/Rep/Ann/"${Ref1Label}".* "${BaseDir}"/Ref/Ann &
#cp "${BaseDir}"/Rep/Asm/WE.genome.fa "${BaseDir}"/Ref/Asm &
#cp "${BaseDir}"/Rep/Asm/CBS432.genome.fa "${BaseDir}"/Ref/Asm &
#cp "${BaseDir}"/Rep/Asm/Jureii.genome.fa "${BaseDir}"/Ref/Asm &

#wait

### part I (run on hulk)

##### WE indexing 
#(

#/usr/bin/time -v bash Index.BWA.sh "${BaseDir}" > Logs/Index.BWA.out 2> Logs/Time.Index.BWA.err &

#/usr/bin/time -v bash Index.SAMt.sh "${BaseDir}" &> Logs/Index.SAMt.log &

#wait
#)

#/usr/bin/time -v bash BWA.Mapping.sh "${BaseDir}" > Logs/BWA.Mapping.out 2> Logs/Time.BWA.Mapping.err

#wait

#/usr/bin/time -v bash SAMt.Dup.sh "${BaseDir}" > Logs/SAMt.Dup.out 2> Logs/Time.SAMt.Dup.err

#wait

#/usr/bin/time -v bash SAMt.Call.sh "${Ref1Label}" "${BaseDir}" > Logs/SAMt.Call.out 2> Logs/Time.SAMt.Call.err &

#wait

#### WE scaffod --> custom genome 

#/usr/bin/time -v Rscript GenGen.R "${BaseDir}" > Logs/GenGen.out 2> Logs/Time.GenGen.err &

#wait

### part II (run on profx)

#### multiple whole-genome alignment WE - samp - SpEU - SJu
#/usr/bin/time -v bash proMAUVEandSNP.sh "${BaseDir}" > Logs/proMAUVEandSNP.out 2> Logs/Time.proMAUVEandSNP.err &
#wait

### Before running the following part rm -r ${BaseDir}"/AdmixtoolsVCF/*"

#rm -r ${BaseDir}"/AdmixtoolsVCF/"

#wait

#Nruns=5
#for j in $(ls -Sr $BaseDir/SNP/*)
#do
# ((i++))
#      (
#       #### create vcf for ADMXTOOLS
#       /usr/bin/time -v Rscript AdToolVCFs.R "${BaseDir}" "${j}"
#      ) &
#if (( i % Nruns == 0 )); then
# wait -n
# i=Nruns-1
#fi
#done

#wait

### replace chrs with arabic numbers 1->16
#/usr/bin/time -v Rscript CHRSplink.R "${BaseDir}" > Logs/CHRSplink.out 2> Logs/Time.CHRSplink.err &

#wait

### put the header

#/usr/bin/time -v bash GenVCF.sh "${BaseDir}" > Logs/GenVCF.out 2> Logs/Time.GenVCF.err &

#wait

#Nruns=5
#for j in $(ls -Sr $BaseDir/AdmixtoolsVCF/*.h.vcf)
#do
# ((i++))
#      (
#       #### gzip samp.h.vcf
#       gzip $j
#       ) &
#if (( i % Nruns == 0 )); then
# wait -n
# i=Nruns-1
#fi
#done

#wait

Nruns=5
for j in $(ls -Sr $BaseDir/AdmixtoolsVCF/*.h.vcf.gz)
do

f=$(echo $j | rev | cut -d"/" -f1 | rev | cut -d"." -f1-2)
f=$(echo $f".snp")

test=$(grep -w $f $BaseDir/done)

if [ -z "$test" ]
  then 
 ((i++))
(
  #### create vcf for ADMXTOOLS
 /usr/bin/time -v bash fromVCF2EIGENSTRAT.sh "${BaseDir}" "$j"
) &
if (( i % Nruns == 0 )); then
 wait -n
 i=Nruns-1
fi
fi
done

#### from vcf to EIGENSTRAT (snp, geno and ind) assuming uniform recomb. rate of 0.3 cM/kb  
#/usr/bin/time -v bash fromVCF2EIGENSTRAT.sh "${BaseDir}" > Logs/fromVCF2EIGENSTRAT.out 2> Logs/Time.fromVCF2EIGENSTRAT.err &

wait

/usr/bin/time -v bash Genind.sh "${BaseDir}" > Logs/Genind.out 2> Logs/Time.Genind.err &

wait

#### performs abba-baba test (Dvalue,Zscore,and sd) for the quartet WE - samp - SpEU - SJu (admixr)
### The abba baba test false-genome vs real genome shows that abba sites are strongly underestimated with the false genome.
### This is due to 1) reduced capacity of performing mapping and calling in corrispondence of introgressions. 2) not optimal capacity of the multiple genomes alignents to confidently detect the introgressed positions.  
### It is raccomanded to use the test below (df) for better understanding of the results.   
/usr/bin/time -v Rscript abbababaT.R "${BaseDir}" > Logs/abbababaT.out 2> Logs/Time.abbababaT.err &

#wait 

#### prepare vcf for estimating f (introgression fraction) and chr-by-chr vcfs for df in 1-kn non-overlapping window.
### about f : f statistic Simon H Martin, Kanchon K Dasmahapatra, Nicola J Nadeau, et al. (2013). Genome-wide evidence for speciation with gene flow in Heliconius butterflies. Genome Res. doi:10.1101/gr.159426.113
### about df : Bastian Pfeifer and Durrell D. Kapan (2019). Estimates of introgression as a function of pairwise distances. BMC Bioinformatics. https://doi.org/10.1186/s12859-019-2747-z
### note : f is reliable if the end-to-end assembly (real genome) of the sample is available.
### note : df overcomes the limitation of f and allows the identification of introgressed genomic windows even if the de novo assembly is not available and it is necessary to create a false genome starting from short reads 

#/usr/bin/time -v bash preFandFmeasure.sh "${BaseDir}" > Logs/preFandFmeasure.out 2> Logs/Time.preFandFmeasure.err &

#wait

#/usr/bin/time -v bash preDFandDFmeasure.sh "${BaseDir}" > Logs/preDFandDFmeasure.out 2> Logs/Time.preDFandDFmeasure.err &
