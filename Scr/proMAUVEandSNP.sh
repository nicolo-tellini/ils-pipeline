Nprocesses=5

BaseDir=$1
GenomeDir=${BaseDir}"/FastaGenomes"
OutDiraln=${BaseDir}"/WGA"
OutDirsnps=${BaseDir}"/SNP"

if [[ ! -d $OutDiraln ]]; then
  mkdir -p $OutDiraln
fi

if [[ ! -d $OutDirsnps ]]; then
  mkdir -p $OutDirsnps
fi

cd ${BaseDir}"/FastaGenomes"

for j in $(ls)
do
if (( i % Nprocesses == 0 )); then
 wait
fi
  ((i++))
(
samp=$(echo $j | cut -d"." -f1)

mkdir -p ${OutDiraln}"/"$samp

cd ${OutDiraln}"/"$samp

control=$(grep -w $samp $BaseDir/already_done.txt)

if [ -z "$control" ]
then

wait
cp -r "${BaseDir}"/Rep/Asm/*fa ${OutDiraln}"/"$samp
cp -r ${GenomeDir}"/"$samp".genome.fa" ${OutDiraln}"/"$samp
wait

aln=WE-"${samp}"-SpEU-SJu
progressiveMauve --output=${OutDiraln}"/"$samp"/"${aln}.xmfa WE.genome.fa $samp".genome.fa" CBS432.genome.fa Jureii.genome.fa

wait

MAUVE_DIR=/scratch/bin/alignment/mauve/
export CLASSPATH="$(find "$MAUVE_DIR" -name \*.jar -print0 | tr '\0' :)$CLASSPATH"
java org.gel.mauve.analysis.SnpExporter -f ${OutDiraln}"/"$samp"/"${aln}.xmfa -o "${OutDirsnps}"/"${aln}".snps

wait

rm ${OutDiraln}"/"$samp"/"*sslist
rm ${OutDiraln}"/"$samp"/"*.genome.fa

fi

)&
done
