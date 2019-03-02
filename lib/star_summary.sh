#!/bin/bash

PREFIX="$1"
SUFFIX="_Log.final.out"
ARRAY=$(ls $1/*"_Log.final.out")
echo -e 'Sample\tReads\tUniquely Mapped%\tUniquely Mapped\tSplices\tMismatch Rate\tDeletion Rate\tInsertion Rate\tMultiple rLoci\tUnmapped\tLink' > $2
echo "Array size: ${#ARRAY[*]}"
echo "Array items:"

for item in $ARRAY
do
    SAMPLE=${item#"$PREFIX"}
    SAMPLE=${SAMPLE%"$SUFFIX"}
    if [[ ${SAMPLE:0:1} == "/" ]]; then SAMPLE=${SAMPLE#"/"}; fi
    OUTPUT1=$(cat $item | grep "Number of input reads" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT2=$(cat $item | grep "Uniquely mapped reads %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT3=$(cat $item | grep "Uniquely mapped reads number" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT4=$(cat $item | grep "Number of splices: Total" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT5=$(cat $item | grep "Mismatch rate per base, %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT6=$(cat $item | grep "Deletion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT7=$(cat $item | grep "Insertion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT8=$(cat $item | grep "% of reads mapped to multiple loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT9=$(cat $item | grep "% of reads mapped to too many loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT10=$(cat $item | grep "% of reads unmapped: too many mismatches" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT11=$(cat $item | grep "% of reads unmapped: too short" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT12=$(cat $item | grep "% of reads unmapped: other" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT13=$(echo $OUTPUT8 + $OUTPUT9 | bc)
    OUTPUT14=$(echo $OUTPUT10 + $OUTPUT11 + $OUTPUT12 | bc)
    OUTPUT15="<a href=\"../results_star/"$SAMPLE"_Log.final.out\" target=\"_blank\">+</a>"
    echo -e ''$SAMPLE'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT12}'%\t'${OUTPUT13}'%\t'$OUTPUT14'%\t'$OUTPUT15 >> $2
done
