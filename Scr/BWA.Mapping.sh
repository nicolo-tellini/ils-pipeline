BaseDir=$1
Nthreads=8

cd $BaseDir

OutDir="MapRaw"

if [[ ! -d $OutDir ]]; then
  mkdir $OutDir
fi

for IndS in $(ls Exp/*gz | cut -d"." -f1 | cut -d"/" -f2 | sort | uniq)
do
  for IndR in ${BaseDir}"/Ref/"*.fa
  do
     if (( i % Nthreads == 0 )); then
       wait
      fi
    ((i++))
  (
   ( RefID=$(basename $IndR | cut -d "." -f 1)
    bwa mem -M -t $Nthreads $IndR Exp/$IndS".R1"*".gz" Exp/$IndS".R2"*".gz" > $OutDir/$IndS.$RefID.sam
    samtools sort --threads $Nthreads -o $OutDir/$IndS.$RefID.srt.bam $OutDir/$IndS.$RefID.sam
    samtools index $OutDir/$IndS.$RefID.srt.bam
    rm -f $OutDir/$IndS.$RefID.sam )&
    ) &
 done
done

wait
