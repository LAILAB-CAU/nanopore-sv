Four methods has been used to detect SVs in B73: Assembly-based method AnchorWave (https://github.com/baoxingsong/anchorwave), ONT-based, Pacbio-based, and short-read sequencing (SRS)-based methods. 

# 1. AnchorWave

 B73_v4 (downloaded from https://maizegdb.org/genome/assembly/Zm-B73-REFERENCE-GRAMENE-4.0) and Mo17_CAU (downloaded from https://maizegdb.org/genome/assembly/Zm-Mo17-REFERENCE-CAU-1.0) were collected. 

Step1. Extract the CDS sequence 

```
anchorwave gff2seq -r ../Zm-Mo17-REFERENCE-CAU-1.0.fa -i ../Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3  -o Zm-Mo17_cds.fa
```

Step2: Align CDS to the ref and query genome respectively
```
minimap2 -x splice -t 6 -k 12 -a -p 0.4 -N 20 ../Zm-B73-REFERENCE-GRAMENE-4.0.fa Zm-Mo17_cds.fa > B73_to_Mo17cds.sam
minimap2 -x splice -t 6 -k 12 -a -p 0.4 -N 20 ../Zm-Mo17-REFERENCE-CAU-1.0.fa Zm-Mo17_cds.fa > Mo17_to_Mo17ref.sam
```

Step3: Genome alignment by anchorWave

```
anchorwave genoAli -i ../Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 -as Zm-Mo17_cds.fa -r ../Zm-Mo17-REFERENCE-CAU-1.0.fa -a B73_to_Mo17cds.sam -ar Mo17_to_Mo17ref.sam -s ../Zm-B73-REFERENCE-GRAMENE-4.0.fa -n B73_Mo17_anchors -o B73_to_Mo17anchorwave.maf -f B73_toMo17anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 12 -IV >B73_Mo17_anchorwave.log 2>&1
```

Step4: Call SV from *.maf file, follow the pipeline described in https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/UserInstructions/CreatePHG_step2_MAFToGVCFPluginDetails.md

```
docker pull maizegenetics/phg
docker pull maizegenetics/phg_liquibase

WORKING_DIR=/data/bxin/nanopore_sv/maf_to_vcf
DOCKER_CONFIG_FILE=/phg/config.txt
cd $WORKING_DIR

echo "MAFToGVCFPlugin.referenceFasta=/phg/mo17.fa" > config.txt
echo "MAFToGVCFPlugin.mafFile=/phg/B73_to_Mo17anchorwave.maf" >>config.txt
echo "MAFToGVCFPlugin.gvcfOutput=/phg/B73_to_Mo17anchorwave.gvcf" >>config.txt 
echo "MAFToGVCFPlugin.sampleName=mo17_anchorwave" >> config.txt
echo "MAFToGVCFPlugin.fillGaps=false" >> config.txt
echo "MAFToGVCFPlugin.twoGvcfs=false" >> config.txt

# The the config file must contain the necessary MAFToGVCFPLugin parameters
docker run --name anchorwave_assemblies_maf --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -MAFToGVCFPlugin -endPlugin

## Seperate DEL and INS
grep -v "#" B73_to_Mo17anchorwave.gvcf | awk 'length($4)>=50||length($5)>=60{print $_}' - > B73_to_Mo17anchorwave.gvcf.temp
python3 convert_gvcf_to_sv.py B73_to_Mo17anchorwave.gvcf.temp  
rm B73_to_Mo17anchorwave.gvcf.temp

## Remove SVs with N in the sequecne 
grep -v N B73_to_Mo17anchorwave.gvcf.del > anchorwave_DEL.bed 
grep -v N B73_to_Mo17anchorwave.gvcf.ins > anchorwave_INS.bed 

## convert bed file to vcf file with SURVIVOR 
SURVIVOR bedtovcf anchorwave_DEL.bed DEL anchorwave_DEL.vcf
SURVIVOR bedtovcf anchorwave_INS.bed INS anchorwave_INS.vcf
```

# 2. ONT-based method

The details in discovering ONT-based SVs has been described in the 'ONT long-read data processing' section. Both Sniffles and PBSV were used to call SVs, and the intersect of Sniffles SVs and PBSV SVs were the final ONT-based SVs, named as intersect_DEL.vcf and intersect_INS.vcf.

# 3. Pacbio-based SVs

The B73 Pacbio data was downloaded from previously-published work at doi:10.1038/nature22971. 

# 4. SRS-based method

We first aligned the SRS reads (B73_1.fq.fastp.gz and B73_2.fq.fastp.gz) by bwa mem and applied both Delly and lumpy in calling SRS-based SVs. 

```
bwa mem -t 60 -R "@RG\tID:B73\tSM:B73\tLB:B73" Zm-Mo17-REFERENCE-CAU-1.0.fa B73_1.fq.fastp.gz B73_2.fq.fastp.gz | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -b - > B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam
```

Run Delly. We used delly call to directly call the SVs with the parameters of -t DEL,INS -g -o, where the -o refers to the genotyped bcf file of SVs.
```
samtools view -b -q 20 B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam | samtools sort -m 3G -@ 10 -o B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam -

samtools index B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam

delly call -t DEL,INS -g Zm-Mo17-REFERENCE-CAU-1.0.fa -o B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.bcf B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam

perl filter_the_InDels_length.pl B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.vcf B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.INS.vcf B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.DEL.vcf
```

Run Lumpy. We extracted the discordant paired-end alignments and split-read alignments and sorted these bam files. Then, we used lumpyexpress to call the raw SVs with the parameters of -B -S -D -O, where the -O refers to the ungenotyped vcf file of SVs. We then used svtyper to genotypes these SVs by the parameters of -i -o and -B. 

```
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam >./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam | extractSplitReads_BwaMem -i stdin | samtools view -b - >./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.splitters.unsorted.bam

# Sort both alignments
samtools sort -@ 16 -o ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.discordants.bam ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.discordants.unsorted.bam

samtools sort -@ 16 -o ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.splitters.bam ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.splitters.unsorted.bam

lumpyexpress -B ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam -S ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.splitters.bam -D ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.discordants.bam -o ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.vcf

samtools sort -@ 10 -m 3G -o ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.bam

samtools index ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam

svtyper -i ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.vcf -o ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.vcf -B ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.bam 

perl ../filter_the_InDels_length.pl ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.vcf ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.INS.vcf ./B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.DEL.vcf
```

Merge Delly and Lumpy results as the final file for SRS-based SVs.

```
echo "B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.DEL.vcf" > merge_delly_lumpy_del.lst
echo "B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.DEL.vcf" >> merge_delly_lumpy_del.lst

echo "B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.sorted.delly.INS.vcf" > merge_delly_lumpy_ins.lst
echo "B73_1.fq.fastp.bwa_mem.Zm-Mo17-REFERENCE-CAU-1.0.lumpy.svtyper.INS.vcf" >> merge_delly_lumpy_ins.lst

SURVIVOR merge merge_delly_lumpy_del.lst 1000 1 1 -1 -1 -1 NGS_DEL.vcf
SURVIVOR merge merge_delly_lumpy_ins.lst 1000 1 1 -1 -1 -1 NGS_INS.vcf
```

# 5. Check the overlap between different methods by SURVIVOR merge

SVs with breakpoints within 1kbp were merged as one SV. 

```
for marker in DEL INS
do
echo "./section2/anchorwave_$marker.vcf" > ./section2/b73_$marker\_123.lst
echo "./section2/intersect_$marker.vcf" >> ./section2/b73_$marker\_123.lst
echo "./section2/NGS_$marker.vcf" >> ./section2/b73_$marker\_123.lst

echo "./section2/anchorwave_$marker.vcf" > ./section2/b73_$marker\_12.lst
echo "./section2/intersect_$marker.vcf" >> ./section2/b73_$marker\_12.lst

echo "./section2/intersect_$marker.vcf" > ./section2/b73_$marker\_23.lst
echo "./section2/NGS_$marker.vcf" >> ./section2/b73_$marker\_23.lst

echo "./section2/anchorwave_$marker.vcf" > ./section2/b73_$marker\_13.lst
echo "./section2/NGS_$marker.vcf" >> ./section2/b73_$marker\_13.lst

# check intersect of methods
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_123.lst 1000 3 1 -1 -1 -1 ./section2/$marker\_123_overlap.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_12.lst 1000 2 1 -1 -1 -1 ./section2/$marker\_12_overlap.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_13.lst 1000 2 1 -1 -1 -1 ./section2/$marker\_13_overlap.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_23.lst 1000 2 1 -1 -1 -1 ./section2/$marker\_23_overlap.vcf

# check the union of methods
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_123.lst 1000 1 1 -1 -1 -1 ./section2/$marker\_123_total.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_12.lst 1000 1 1 -1 -1 -1 ./section2/$marker\_12_total.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_13.lst 1000 1 1 -1 -1 -1 ./section2/$marker\_13_total.vcf
/public1/home/sc30797/bxin/SURVIVOR/Debug/SURVIVOR merge ./section2/b73_$marker\_23.lst 1000 1 1 -1 -1 -1 ./section2/$marker\_23_total.vcf

echo $marker
assembly_NGS_ONT=`grep -v "#" ./section2/$marker\_123_overlap.vcf | wc -l`
assembly_NGS=`grep -v "#" ./section2/$marker\_13_overlap.vcf | wc -l`
assembly_ONT=`grep -v "#" ./section2/$marker\_12_overlap.vcf | wc -l`
NGS_ONT=`grep -v "#" ./section2/$marker\_23_overlap.vcf | wc -l`
echo "$assembly_ONT $NGS_ONT $assembly_NGS  $assembly_NGS_ONT"
grep -v "#" ./section2/$marker\_123_total.vcf | wc -l
grep -v "#" ./section2/$marker\_12_total.vcf | wc -l
grep -v "#" ./section2/$marker\_13_total.vcf | wc -l
grep -v "#" ./section2/$marker\_23_total.vcf | wc -l
done
```