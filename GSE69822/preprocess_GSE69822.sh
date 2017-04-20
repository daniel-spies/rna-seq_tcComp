outDir=/cclab_nas/dspies/data/simulated/study/data/GSE69822
toolDir=~/software/sratoolkit.2.6.3-ubuntu64/bin
gtf=/cclab_nas/Reference/ensembl/Homo_sapiens.GRCh38.83.gtf
genDir=/cclab_nas/Reference/star/GRCh38_83

## downloading data
for file in $(cat $outDir/SRR_Acc_List.txt)
do
    $toolDir/fastq-dump -O $outDir/raw $file &
done

# mapping data ~ 89% uniquely mapped reads
for file in $(ls $outDir/raw)
do
    STAR --readFilesIn $outDir/raw/$file  --outFileNamePrefix $outDir/mapped/${file}_ \
        --outFilterMatchNmin 20 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.02 --outFilterMultimapNmax 3000 \
        --genomeDir $genDir --sjdbGTFfile $gtf --outFilterIntronMotifs RemoveNoncanonical --outSAMmode Full \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 40000000000 \
        --runThreadN 25 --outStd Log
done

# quantify reads
featureCounts -T 30 -a $gtf -s 0 -o $outDir/GSE69822.counts $(ls $outDir/mapped/*Aligned.sortedByCoord.out.bam)